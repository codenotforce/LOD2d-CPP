#pragma once

#include "helmholtz/model.h"
#include "helmholtz/schwarz_patch.h"
#include "helmholtz/schwarz_local_solver.h"

#include <memory>
#include <vector>

namespace lod2d::helmholtz {

enum class HelmholtzTwoLevelSchwarzMode {
    Additive,
    Hybrid
};

enum class HelmholtzSchwarzExtension {
    WeightedOverlap,
    RestrictedCore
};

enum class HelmholtzSchwarzFactorizationReuse {
    None,
    IdenticalMatrix
};

struct HelmholtzTwoLevelSchwarzConfig {
    HelmholtzSchwarzArtificialBoundary artificial_boundary =
        HelmholtzSchwarzArtificialBoundary::HomogeneousDirichlet;
    double artificial_impedance_beta = 1.0;
    HelmholtzSchwarzExtension extension =
        HelmholtzSchwarzExtension::WeightedOverlap;
    HelmholtzSchwarzFactorizationReuse factorization_reuse =
        HelmholtzSchwarzFactorizationReuse::None;
    HelmholtzSchwarzLocalSolverConfig local_solver;
};

struct HelmholtzSchwarzLocalSolverDiagnostics {
    long long solve_calls = 0;
    long long total_iterations = 0;
    long long total_restarts = 0;
    int max_iterations = 0;
    double max_relative_residual = 0.0;
    int max_vcycle_levels = 0;
    int min_vcycle_coarse_dofs = 0;
    int max_vcycle_finest_dofs = 0;
};

struct HelmholtzTwoLevelSchwarzDiagnostics {
    int subdomains = 0;
    int min_local_dofs = 0;
    int max_local_dofs = 0;
    int min_owned_dofs = 0;
    int max_owned_dofs = 0;
    int uncovered_dofs = 0;
    double partition_unity_error = 0.0;
    int artificial_boundary_edges = 0;
    int physical_boundary_edges = 0;
    int local_solver_groups = 0;
    int reused_factorizations = 0;
    int max_reuse_group = 0;
};

class HelmholtzTwoLevelSchwarzPreconditioner {
public:
    explicit HelmholtzTwoLevelSchwarzPreconditioner(
        const HelmholtzLodModel &model,
        HelmholtzTwoLevelSchwarzConfig config = {});
    HelmholtzTwoLevelSchwarzPreconditioner(HelmholtzLodModel &&) = delete;
    HelmholtzTwoLevelSchwarzPreconditioner(
        HelmholtzLodModel &&,
        HelmholtzTwoLevelSchwarzConfig) = delete;

    ComplexVector apply_coarse(const ComplexVector &right_hand_side) const;
    ComplexVector apply_local(const ComplexVector &right_hand_side) const;
    ComplexVector apply(
        const ComplexVector &right_hand_side,
        HelmholtzTwoLevelSchwarzMode mode) const;

    const HelmholtzTwoLevelSchwarzDiagnostics &diagnostics() const {
        return diagnostics_;
    }
    const HelmholtzTwoLevelSchwarzConfig &config() const { return config_; }
    HelmholtzSchwarzLocalSolverDiagnostics local_solver_diagnostics() const {
        return local_solver_diagnostics_;
    }

private:
    struct Subdomain {
        std::vector<int> global_dofs;
        Eigen::VectorXd injection_weights;
    };

    struct SolverGroup {
        std::vector<int> subdomains;
        std::unique_ptr<HelmholtzSchwarzLocalSolver> solver;
    };

    void validate_right_hand_side(const ComplexVector &right_hand_side) const;

    HelmholtzTwoLevelSchwarzConfig config_;
    const HelmholtzLodModel *model_ = nullptr;
    const ComplexSparseMatrix *fine_operator_ = nullptr;
    int dimension_ = 0;
    std::vector<Subdomain> subdomains_;
    std::vector<SolverGroup> solver_groups_;
    HelmholtzTwoLevelSchwarzDiagnostics diagnostics_;
    mutable HelmholtzSchwarzLocalSolverDiagnostics local_solver_diagnostics_;
};

} // namespace lod2d::helmholtz
