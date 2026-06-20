#pragma once
#include "lod/corrector.h"
#include "mesh/types.h"
#include "solver/lod_reuse.h"
#include <Eigen/Sparse>
#include <array>
#include <memory>
#include <vector>

namespace lod2d {

struct LodProblemConfig {
    int H = 4;
    int h = 10;
    int ell = 2;
    int d = 2;
    CorrectorSolver solver = CorrectorSolver::Cholmod;
    bool keep_setup_matrices = false;
    TriMesh initial_mesh;
};

struct LodProblemData {
    TriMesh coarse;
    TriMesh fine;
    int NH = 0;
    int NTH = 0;
    int Nh = 0;
    int NTh = 0;
    std::vector<int> nngH;
    std::vector<int> nngh;
    std::vector<std::array<int, 3>> dghidx;
    std::vector<std::array<int, 3>> dgHidx;
    Eigen::SparseMatrix<double> P_node;
    Eigen::SparseMatrix<double> P_elem;
    Eigen::SparseMatrix<double> P_dg;
};

struct LodOperators {
    ElementStiffnessBlocks element_stiffness;
    Eigen::SparseMatrix<double> Sh;
    Eigen::SparseMatrix<double> Mh;
    Eigen::SparseMatrix<double> patch;
    InterpolationRows interpolation_rows;
    FineElementChildren fine_element_children;
};

TriMesh make_unit_square_mesh();

LodProblemData build_lod_problem_data(const TriMesh &initial_mesh, int H, int h);

LodOperators build_lod_operators(
    const LodProblemData &problem,
    const std::vector<double> &Ah,
    int ell);

std::vector<CorrectorEntries> build_lod_correctors(
    const LodProblemData &problem,
    const LodOperators &operators,
    int d,
    CorrectorSolver solver);

Eigen::SparseMatrix<double> build_lod_basis(
    const LodProblemData &problem,
    const std::vector<CorrectorEntries> &correctors);

class LodModel {
public:
    LodModel() = default;

    static LodModel build(const LodProblemConfig &config, const std::vector<double> &Ah);

    const LodProblemConfig &config() const { return config_; }
    const LodProblemData &problem() const { return problem_; }
    const LodReusableSystem &reusable_system() const;

    LodReuseSolution solve_from_coarse_values(const Eigen::VectorXd &f_coarse) const;
    LodReuseSolution solve_from_fine_values(const Eigen::VectorXd &f_fine) const;

private:
    LodProblemConfig config_;
    LodProblemData problem_;
    std::unique_ptr<LodReusableSystem> system_;
};

} // namespace lod2d
