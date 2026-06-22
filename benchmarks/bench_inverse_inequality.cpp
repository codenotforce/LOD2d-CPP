/// Numerical check for the LOD inverse inequality on the multiscale space.
#include "lod/lod_model.h"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lod2d;
namespace chr = std::chrono;

namespace {

struct Options {
    int H = 4;
    int h = 9;
    int ell = 2;
    int d = 2;
    int threads = -1;
    int h_minus_H = 5;
    int H_min = 2;
    int H_max = 4;
    int h_min = 0;
    int h_max = 0;
    int ell_min = 0;
    int ell_max = 0;
    bool sweep_H = false;
    bool sweep_h = false;
    bool sweep_ell = false;
    bool free_space = true;
    CorrectorSolver solver = CorrectorSolver::EigenLLT;
    std::string coeff = "unit";
    std::string solver_spec = "auto";
    std::string basis = "lod";
    std::string numerator = "lod";
    std::string cache_dir = "results/corrector_cache";
};

struct Stats {
    double min = 0.0;
    double median = 0.0;
    double p90 = 0.0;
    double p99 = 0.0;
    double max = 0.0;
    int argmax = -1;
};

const char *solver_name(CorrectorSolver solver) {
    if (solver == CorrectorSolver::Cholmod) return "cholmod";
    if (solver == CorrectorSolver::CholmodCached) return "cholmod_cached";
    return "eigen";
}

CorrectorSolver parse_solver(const std::string &value, int h) {
    if (value == "eigen") return CorrectorSolver::EigenLLT;
    if (value == "cholmod") return CorrectorSolver::Cholmod;
    if (value == "cholmod_cached") return CorrectorSolver::CholmodCached;
    if (value == "auto") return (h >= 10) ? CorrectorSolver::Cholmod : CorrectorSolver::EigenLLT;
    throw std::invalid_argument("unknown solver: " + value);
}

Options parse_options(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--H=", 0) == 0) opt.H = std::stoi(arg.substr(4));
        else if (arg.rfind("--h=", 0) == 0) opt.h = std::stoi(arg.substr(4));
        else if (arg.rfind("--ell=", 0) == 0) opt.ell = std::stoi(arg.substr(6));
        else if (arg.rfind("--coeff=", 0) == 0) opt.coeff = arg.substr(8);
        else if (arg.rfind("--solver=", 0) == 0) opt.solver_spec = arg.substr(9);
        else if (arg.rfind("--basis=", 0) == 0) opt.basis = arg.substr(8);
        else if (arg.rfind("--numerator=", 0) == 0) opt.numerator = arg.substr(12);
        else if (arg.rfind("--cache-dir=", 0) == 0) opt.cache_dir = arg.substr(12);
        else if (arg == "--no-cache") opt.cache_dir.clear();
        else if (arg.rfind("--threads=", 0) == 0) {
            std::string value = arg.substr(10);
            if (value == "auto") opt.threads = -1;
            else if (value == "env") opt.threads = 0;
            else opt.threads = std::stoi(value);
        } else if (arg == "--sweep-H") opt.sweep_H = true;
        else if (arg == "--sweep-h") opt.sweep_h = true;
        else if (arg == "--sweep-ell") opt.sweep_ell = true;
        else if (arg.rfind("--h-minus-H=", 0) == 0) opt.h_minus_H = std::stoi(arg.substr(12));
        else if (arg.rfind("--H-min=", 0) == 0) opt.H_min = std::stoi(arg.substr(8));
        else if (arg.rfind("--H-max=", 0) == 0) opt.H_max = std::stoi(arg.substr(8));
        else if (arg.rfind("--h-min=", 0) == 0) opt.h_min = std::stoi(arg.substr(8));
        else if (arg.rfind("--h-max=", 0) == 0) opt.h_max = std::stoi(arg.substr(8));
        else if (arg.rfind("--ell-min=", 0) == 0) opt.ell_min = std::stoi(arg.substr(10));
        else if (arg.rfind("--ell-max=", 0) == 0) opt.ell_max = std::stoi(arg.substr(10));
        else if (arg == "--space=free") opt.free_space = true;
        else if (arg == "--space=all") opt.free_space = false;
        else {
            throw std::invalid_argument(
                "usage: bench_inverse_inequality [--H=N --h=N --ell=N] "
                "[--solver=eigen|cholmod|cholmod_cached|auto] "
                "[--coeff=unit|file:PATH|checkerboard:CONTRAST] "
                "[--basis=lod|coarse] [--numerator=lod|corrector] [--cache-dir=PATH|--no-cache] [--space=free|all] [--threads=auto|env|N] "
                "[--sweep-H --H-min=N --H-max=N --h-minus-H=N] "
                "[--sweep-h --h-min=N --h-max=N] "
                "[--sweep-ell --ell-min=N --ell-max=N]");
        }
    }
    if (opt.H < 0 || opt.h < opt.H) throw std::invalid_argument("require 0 <= H <= h");
    if (opt.ell < 0) throw std::invalid_argument("ell must be nonnegative");
    if (opt.sweep_H && opt.H_max < opt.H_min) throw std::invalid_argument("require H-min <= H-max");
    const int sweep_count = (opt.sweep_H ? 1 : 0) + (opt.sweep_h ? 1 : 0) + (opt.sweep_ell ? 1 : 0);
    if (sweep_count > 1) throw std::invalid_argument("choose only one sweep mode");
    if (opt.sweep_h) {
        if (opt.h_min == 0) opt.h_min = opt.H + 1;
        if (opt.h_max == 0) opt.h_max = opt.h;
        if (opt.h_min < opt.H || opt.h_max < opt.h_min)
            throw std::invalid_argument("require H <= h-min <= h-max");
    }
    if (opt.sweep_ell) {
        if (opt.ell_min == 0) opt.ell_min = 1;
        if (opt.ell_max == 0) opt.ell_max = opt.ell;
        if (opt.ell_min < 0 || opt.ell_max < opt.ell_min)
            throw std::invalid_argument("require 0 <= ell-min <= ell-max");
    }
    if (opt.basis != "lod" && opt.basis != "coarse") throw std::invalid_argument("basis must be lod or coarse");
    if (opt.numerator != "lod" && opt.numerator != "corrector") throw std::invalid_argument("numerator must be lod or corrector");
    opt.solver = parse_solver(opt.solver_spec, opt.h);
    return opt;
}

#ifdef _OPENMP
void apply_thread_option(const Options &opt, int h) {
    if (opt.threads > 0) {
        omp_set_num_threads(opt.threads);
        return;
    }
    if (opt.threads == 0) return;
    if (h >= 10 && (opt.solver == CorrectorSolver::Cholmod || opt.solver == CorrectorSolver::CholmodCached) && omp_get_max_threads() > 8)
        omp_set_num_threads(8);
}
#else
void apply_thread_option(const Options &, int) {}
#endif

std::vector<double> read_coeff_file(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open coefficient file: " + path);
    int n = 0;
    in >> n;
    std::vector<double> coeff(n);
    for (int i = 0; i < n; ++i) in >> coeff[i];
    if (!in) throw std::runtime_error("failed to read coefficient file: " + path);
    return coeff;
}

std::vector<double> make_coefficients(const LodProblemData &problem, const std::string &spec) {
    if (spec == "unit") return std::vector<double>(problem.NTh, 1.0);
    if (spec.rfind("file:", 0) == 0) {
        std::vector<double> coeff = read_coeff_file(spec.substr(5));
        if (static_cast<int>(coeff.size()) != problem.NTh)
            throw std::runtime_error("coefficient file length does not match fine element count");
        return coeff;
    }
    if (spec.rfind("checkerboard:", 0) == 0) {
        const double contrast = std::stod(spec.substr(13));
        std::vector<double> coeff(problem.NTh, 1.0);
        const int cells = 8;
        for (int e = 0; e < problem.NTh; ++e) {
            Eigen::Vector2d c = Eigen::Vector2d::Zero();
            for (int v : problem.fine.elems[e]) c += problem.fine.nodes[v];
            c /= 3.0;
            const int ix = std::min(cells - 1, std::max(0, static_cast<int>(std::floor(cells * c.x()))));
            const int iy = std::min(cells - 1, std::max(0, static_cast<int>(std::floor(cells * c.y()))));
            coeff[e] = ((ix + iy) % 2 == 0) ? contrast : 1.0;
        }
        return coeff;
    }
    throw std::invalid_argument("unknown coefficient spec: " + spec);
}


std::string sanitize_key(std::string value) {
    for (char &c : value) {
        const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                        (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.';
        if (!ok) c = '_';
    }
    return value;
}

std::string cache_path(const Options &opt, const LodProblemData &problem) {
    std::ostringstream name;
    name << "ct_H" << opt.H
         << "_h" << opt.h
         << "_ell" << opt.ell
         << "_d" << opt.d
         << "_Nh" << problem.Nh
         << "_NTh" << problem.NTh
         << "_coeff_" << sanitize_key(opt.coeff)
         << "_solver_" << solver_name(opt.solver)
         << ".bin";
    return (std::filesystem::path(opt.cache_dir) / name.str()).string();
}

template <class T>
void write_binary(std::ostream &out, const T &value) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <class T>
bool read_binary(std::istream &in, T &value) {
    return static_cast<bool>(in.read(reinterpret_cast<char *>(&value), sizeof(T)));
}

void write_string(std::ostream &out, const std::string &value) {
    const std::uint64_t n = static_cast<std::uint64_t>(value.size());
    write_binary(out, n);
    out.write(value.data(), static_cast<std::streamsize>(n));
}

bool read_string(std::istream &in, std::string &value) {
    std::uint64_t n = 0;
    if (!read_binary(in, n)) return false;
    value.resize(static_cast<size_t>(n));
    return static_cast<bool>(in.read(value.data(), static_cast<std::streamsize>(n)));
}

bool load_corrector_cache(
    const std::string &path,
    const Options &opt,
    const LodProblemData &problem,
    std::vector<CorrectorEntries> &correctors) {
    if (path.empty() || !std::filesystem::exists(path)) return false;
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    std::string magic;
    if (!read_string(in, magic) || magic != "LOD2D_CORRECTOR_CACHE_V1") return false;
    std::string coeff;
    if (!read_string(in, coeff) || coeff != opt.coeff) return false;
    std::string solver;
    if (!read_string(in, solver) || solver != solver_name(opt.solver)) return false;

    int H = 0, h = 0, ell = 0, d = 0, NH = 0, Nh = 0, NTH = 0, NTh = 0;
    if (!read_binary(in, H) || !read_binary(in, h) || !read_binary(in, ell) || !read_binary(in, d) ||
        !read_binary(in, NH) || !read_binary(in, Nh) || !read_binary(in, NTH) || !read_binary(in, NTh))
        return false;
    if (H != opt.H || h != opt.h || ell != opt.ell || d != opt.d ||
        NH != problem.NH || Nh != problem.Nh || NTH != problem.NTH || NTh != problem.NTh)
        return false;

    std::uint64_t n_correctors = 0;
    if (!read_binary(in, n_correctors)) return false;
    correctors.assign(static_cast<size_t>(n_correctors), {});
    for (auto &entries : correctors) {
        std::uint64_t n_entries = 0;
        if (!read_binary(in, n_entries)) return false;
        entries.resize(static_cast<size_t>(n_entries));
        for (auto &entry : entries) {
            if (!read_binary(in, entry.row) || !read_binary(in, entry.col) || !read_binary(in, entry.value))
                return false;
        }
    }
    return static_cast<bool>(in);
}

void save_corrector_cache(
    const std::string &path,
    const Options &opt,
    const LodProblemData &problem,
    const std::vector<CorrectorEntries> &correctors) {
    if (path.empty()) return;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    const std::string tmp_path = path + ".tmp";
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("failed to open corrector cache for write: " + tmp_path);

    write_string(out, "LOD2D_CORRECTOR_CACHE_V1");
    write_string(out, opt.coeff);
    write_string(out, solver_name(opt.solver));
    write_binary(out, opt.H);
    write_binary(out, opt.h);
    write_binary(out, opt.ell);
    write_binary(out, opt.d);
    write_binary(out, problem.NH);
    write_binary(out, problem.Nh);
    write_binary(out, problem.NTH);
    write_binary(out, problem.NTh);
    const std::uint64_t n_correctors = static_cast<std::uint64_t>(correctors.size());
    write_binary(out, n_correctors);
    for (const auto &entries : correctors) {
        const std::uint64_t n_entries = static_cast<std::uint64_t>(entries.size());
        write_binary(out, n_entries);
        for (const auto &entry : entries) {
            write_binary(out, entry.row);
            write_binary(out, entry.col);
            write_binary(out, entry.value);
        }
    }
    if (!out) throw std::runtime_error("failed while writing corrector cache: " + tmp_path);
    out.close();
    std::filesystem::rename(tmp_path, path);
}

Eigen::Matrix3d unit_stiffness(const TriMesh &mesh, int e) {
    const int a = mesh.elems[e][0];
    const int b = mesh.elems[e][1];
    const int c = mesh.elems[e][2];
    const auto &pa = mesh.nodes[a];
    const auto &pb = mesh.nodes[b];
    const auto &pc = mesh.nodes[c];

    const double ve1_x = pc.x() - pb.x(), ve1_y = pc.y() - pb.y();
    const double ve2_x = pa.x() - pc.x(), ve2_y = pa.y() - pc.y();
    const double ve3_x = pb.x() - pa.x(), ve3_y = pb.y() - pa.y();
    const double area2 = std::abs(ve3_x * ve2_y - ve3_y * ve2_x);
    const double area = 0.5 * area2;
    const double inv_a2 = (area2 > 0.0) ? 1.0 / area2 : 0.0;

    const double grad[3][2] = {
        {-ve1_y * inv_a2, ve1_x * inv_a2},
        {-ve2_y * inv_a2, ve2_x * inv_a2},
        {-ve3_y * inv_a2, ve3_x * inv_a2}
    };
    Eigen::Matrix3d S;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            S(i, j) = area * (grad[i][0] * grad[j][0] + grad[i][1] * grad[j][1]);
    return S;
}

Eigen::Matrix3d mass_matrix(const TriMesh &mesh, int e) {
    const auto &pa = mesh.nodes[mesh.elems[e][0]];
    const auto &pb = mesh.nodes[mesh.elems[e][1]];
    const auto &pc = mesh.nodes[mesh.elems[e][2]];
    const double area = 0.5 * std::abs((pb.x() - pa.x()) * (pc.y() - pa.y()) -
                                       (pb.y() - pa.y()) * (pc.x() - pa.x()));
    Eigen::Matrix3d M;
    M << 2.0, 1.0, 1.0,
         1.0, 2.0, 1.0,
         1.0, 1.0, 2.0;
    return (area / 12.0) * M;
}

double coarse_diameter(const TriMesh &mesh, int e) {
    double diameter = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            diameter = std::max(diameter, (mesh.nodes[mesh.elems[e][i]] - mesh.nodes[mesh.elems[e][j]]).norm());
        }
    }
    return diameter;
}

Stats summarize(const std::vector<double> &values) {
    if (values.empty()) return {};
    std::vector<int> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return values[a] < values[b]; });
    auto pick = [&](double q) {
        const int idx = std::min(static_cast<int>(order.size()) - 1,
                                 std::max(0, static_cast<int>(std::ceil(q * order.size()) - 1)));
        return values[order[idx]];
    };
    Stats out;
    out.min = values[order.front()];
    out.median = pick(0.50);
    out.p90 = pick(0.90);
    out.p99 = pick(0.99);
    out.max = values[order.back()];
    out.argmax = order.back();
    return out;
}

struct RunResult {
    int H = 0;
    int h = 0;
    int ell = 0;
    int coarse_elements = 0;
    int fine_elements = 0;
    Stats stats;
    double setup_ms = 0.0;
    double inverse_ms = 0.0;
};

RunResult run_case(Options opt) {
    apply_thread_option(opt, opt.h);
    auto t0 = chr::high_resolution_clock::now();

    LodProblemData problem = build_lod_problem_data(make_unit_square_mesh(), opt.H, opt.h);
    FineElementChildren fine_element_children;
    Eigen::SparseMatrix<double> G;
    if (opt.basis == "coarse") {
        fine_element_children = build_fine_element_children(problem.P_elem, problem.NTH);
        G = problem.P_node;
    } else {
        std::vector<CorrectorEntries> correctors;
        const std::string ct_cache_path = opt.cache_dir.empty() ? std::string() : cache_path(opt, problem);
        const bool cache_hit = load_corrector_cache(ct_cache_path, opt, problem, correctors);
        if (cache_hit) {
            std::cout << "Corrector cache: hit " << ct_cache_path << std::endl;
            fine_element_children = build_fine_element_children(problem.P_elem, problem.NTH);
        } else {
            if (!ct_cache_path.empty()) std::cout << "Corrector cache: miss " << ct_cache_path << std::endl;
            std::vector<double> coeff = make_coefficients(problem, opt.coeff);
            LodOperators operators = build_lod_operators(problem, coeff, opt.ell);
            correctors = build_lod_correctors(problem, operators, opt.d, opt.solver);
            fine_element_children = std::move(operators.fine_element_children);
            save_corrector_cache(ct_cache_path, opt, problem, correctors);
        }
        G = build_lod_basis(problem, correctors);
    }
    auto t1 = chr::high_resolution_clock::now();

    std::vector<char> is_dirichlet(problem.NH, false);
    for (int v : problem.coarse.dirichlet) {
        if (v >= 0 && v < problem.NH) is_dirichlet[v] = true;
    }

    Eigen::SparseMatrix<double> numerator_matrix;
    if (opt.numerator == "corrector") numerator_matrix = problem.P_node - G;
    else numerator_matrix = G;

    auto build_rows = [&](const Eigen::SparseMatrix<double> &matrix) {
        std::vector<std::vector<std::pair<int, double>>> rows(problem.Nh);
        for (int col = 0; col < matrix.outerSize(); ++col) {
            if (opt.free_space && is_dirichlet[col]) continue;
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, col); it; ++it) {
                if (it.value() != 0.0) rows[it.row()].push_back({col, it.value()});
            }
        }
        return rows;
    };

    auto denominator_rows = build_rows(G);
    auto numerator_rows = build_rows(numerator_matrix);

    std::vector<int> fine_stamp(problem.Nh, 0);
    std::vector<int> col_stamp(problem.NH, 0);
    std::vector<int> col_to_local(problem.NH, -1);
    int stamp = 0;
    std::vector<double> ratios(problem.NTH, 0.0);

    for (int T = 0; T < problem.NTH; ++T) {
        ++stamp;
        std::vector<int> local_fine_vertices;
        std::vector<int> active_cols;
        for (int fe : fine_element_children[T]) {
            for (int fv : problem.fine.elems[fe]) {
                if (fine_stamp[fv] != stamp) {
                    fine_stamp[fv] = stamp;
                    local_fine_vertices.push_back(fv);
                    auto touch_col = [&](int col) {
                        if (col_stamp[col] != stamp) {
                            col_stamp[col] = stamp;
                            col_to_local[col] = static_cast<int>(active_cols.size());
                            active_cols.push_back(col);
                        }
                    };
                    for (const auto &[col, value] : denominator_rows[fv]) {
                        (void)value;
                        touch_col(col);
                    }
                    for (const auto &[col, value] : numerator_rows[fv]) {
                        (void)value;
                        touch_col(col);
                    }
                }
            }
        }

        const int n = static_cast<int>(active_cols.size());
        if (n == 0) continue;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);

        for (int fe : fine_element_children[T]) {
            Eigen::Matrix3d Se = unit_stiffness(problem.fine, fe);
            Eigen::Matrix3d Me = mass_matrix(problem.fine, fe);
            Eigen::MatrixXd B_num = Eigen::MatrixXd::Zero(3, n);
            Eigen::MatrixXd B_den = Eigen::MatrixXd::Zero(3, n);
            for (int i = 0; i < 3; ++i) {
                const int fv = problem.fine.elems[fe][i];
                for (const auto &[col, value] : numerator_rows[fv]) {
                    if (col_stamp[col] == stamp) B_num(i, col_to_local[col]) = value;
                }
                for (const auto &[col, value] : denominator_rows[fv]) {
                    if (col_stamp[col] == stamp) B_den(i, col_to_local[col]) = value;
                }
            }
            A.noalias() += B_num.transpose() * Se * B_num;
            M.noalias() += B_den.transpose() * Me * B_den;
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> m_eig(M);
        if (m_eig.info() != Eigen::Success)
            throw std::runtime_error("local mass eigensolve failed");
        const Eigen::VectorXd mass_evals = m_eig.eigenvalues();
        const double mass_max = mass_evals.size() ? mass_evals.maxCoeff() : 0.0;
        if (mass_max <= 0.0) continue;
        const double tol = std::max(1e-30, 1e-12 * mass_max);
        int keep = 0;
        for (int i = 0; i < mass_evals.size(); ++i) {
            if (mass_evals(i) > tol) ++keep;
        }
        if (keep == 0) continue;
        Eigen::MatrixXd Z(n, keep);
        int k = 0;
        for (int i = 0; i < mass_evals.size(); ++i) {
            if (mass_evals(i) > tol) {
                Z.col(k++) = m_eig.eigenvectors().col(i) / std::sqrt(mass_evals(i));
            }
        }
        Eigen::MatrixXd reduced = Z.transpose() * A * Z;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> a_eig(reduced);
        if (a_eig.info() != Eigen::Success)
            throw std::runtime_error("local stiffness eigensolve failed");
        const double lambda_max = std::max(0.0, a_eig.eigenvalues().maxCoeff());
        ratios[T] = coarse_diameter(problem.coarse, T) * std::sqrt(lambda_max);
    }

    auto t2 = chr::high_resolution_clock::now();
    RunResult out;
    out.H = opt.H;
    out.h = opt.h;
    out.ell = opt.ell;
    out.coarse_elements = problem.NTH;
    out.fine_elements = problem.NTh;
    out.stats = summarize(ratios);
    out.setup_ms = chr::duration<double, std::milli>(t1 - t0).count();
    out.inverse_ms = chr::duration<double, std::milli>(t2 - t1).count();
    return out;
}

void print_result(const Options &opt, const RunResult &r) {
    std::cout << "=== Inverse inequality check ===\n";
    std::cout << "H=" << r.H << " h=" << r.h << " ell=" << r.ell
              << " coeff=" << opt.coeff << " solver=" << solver_name(opt.solver)
              << " basis=" << opt.basis
              << " numerator=" << opt.numerator
              << " space=" << (opt.free_space ? "free" : "all") << "\n";
#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n";
#endif
    std::cout << "Coarse elements: " << r.coarse_elements << "  Fine elements: " << r.fine_elements << "\n";
    std::cout << "Setup: " << r.setup_ms << " ms\n";
    std::cout << "Inverse scan: " << r.inverse_ms << " ms\n\n";
    if (opt.numerator == "corrector")
        std::cout << "Q_T = H_T * ||grad C v||_T / ||(1-C)v||_T\n";
    else
        std::cout << "Q_T = H_T * ||grad (1-C)v||_T / ||(1-C)v||_T\n";
    std::cout << "min     : " << r.stats.min << "\n";
    std::cout << "median  : " << r.stats.median << "\n";
    std::cout << "p90     : " << r.stats.p90 << "\n";
    std::cout << "p99     : " << r.stats.p99 << "\n";
    std::cout << "max     : " << r.stats.max << "\n";
    std::cout << "argmax T: " << r.stats.argmax << "\n";
}

void print_sweep_header() {
    std::cout << "H,h,ell,basis,numerator,coarse_elements,fine_elements,min,median,p90,p99,max,argmax,setup_ms,inverse_ms\n";
}

void print_sweep_row(const Options &opt, const RunResult &r) {
    std::cout << r.H << ',' << r.h << ',' << r.ell << ','
              << opt.basis << ',' << opt.numerator << ',' << r.coarse_elements << ',' << r.fine_elements << ','
              << r.stats.min << ',' << r.stats.median << ',' << r.stats.p90 << ','
              << r.stats.p99 << ',' << r.stats.max << ',' << r.stats.argmax << ','
              << r.setup_ms << ',' << r.inverse_ms << "\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        Options opt = parse_options(argc, argv);
        if (!opt.sweep_H && !opt.sweep_h && !opt.sweep_ell) {
            RunResult result = run_case(opt);
            print_result(opt, result);
            return 0;
        }

        print_sweep_header();
        if (opt.sweep_h) {
            for (int h = opt.h_min; h <= opt.h_max; ++h) {
                Options run = opt;
                run.h = h;
                run.solver = parse_solver(opt.solver_spec, run.h);
                RunResult result = run_case(run);
                print_sweep_row(run, result);
            }
            return 0;
        }

        if (opt.sweep_ell) {
            for (int ell = opt.ell_min; ell <= opt.ell_max; ++ell) {
                Options run = opt;
                run.ell = ell;
                RunResult result = run_case(run);
                print_sweep_row(run, result);
            }
            return 0;
        }

        for (int H = opt.H_min; H <= opt.H_max; ++H) {
            Options run = opt;
            run.H = H;
            run.h = H + opt.h_minus_H;
            run.solver = parse_solver(opt.solver_spec, run.h);
            RunResult result = run_case(run);
            print_sweep_row(run, result);
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
