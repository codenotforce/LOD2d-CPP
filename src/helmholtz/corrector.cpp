#include "helmholtz/corrector.h"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lod2d::helmholtz {
namespace {

struct ColumnEntry {
    int row = -1;
    Complex value = 0.0;
};

using ColumnEntries = std::vector<ColumnEntry>;

std::vector<std::vector<std::pair<int, int>>> build_node_incidence(
    const TriMesh &coarse) {
    std::vector<std::vector<std::pair<int, int>>> incidence(
        coarse.nodes.size());
    for (int element = 0;
         element < static_cast<int>(coarse.elems.size());
         ++element) {
        for (int local = 0; local < 3; ++local)
            incidence[coarse.elems[element][local]].emplace_back(element, local);
    }
    return incidence;
}

ColumnEntries assemble_column(
    int column,
    int fine_node_count,
    const std::vector<std::vector<std::pair<int, int>>> &incidence,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    ColumnEntries raw;
    std::size_t upper_bound = 0;
    for (const auto &[element, local] : incidence[column]) {
        (void)local;
        upper_bound += correctors[element].size() / 3 + 1;
    }
    raw.reserve(upper_bound);
    for (const auto &[element, local] : incidence[column]) {
        for (const HelmholtzCorrectorEntry &entry : correctors[element]) {
            if (entry.local_coarse_vertex != local) continue;
            if (entry.row < 0 || entry.row >= fine_node_count)
                throw std::out_of_range(
                    "Helmholtz corrector row lies outside the fine mesh");
            raw.push_back({entry.row, entry.value});
        }
    }
    std::stable_sort(raw.begin(), raw.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.row < rhs.row;
    });

    ColumnEntries merged;
    for (std::size_t first = 0; first < raw.size();) {
        std::size_t last = first + 1;
        Complex value = raw[first].value;
        while (last < raw.size() && raw[last].row == raw[first].row) {
            value += raw[last].value;
            ++last;
        }
        if (value != Complex(0.0, 0.0))
            merged.push_back({raw[first].row, value});
        first = last;
    }
    merged.shrink_to_fit();
    return merged;
}

} // namespace

ComplexSparseMatrix build_helmholtz_corrector_matrix(
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    if (correctors.size() != coarse.elems.size())
        throw std::invalid_argument("Helmholtz corrector count must match coarse elements");
    if (fine_node_count < 0)
        throw std::invalid_argument("Helmholtz fine node count must be nonnegative");

    for (const auto &element : correctors) {
        for (const HelmholtzCorrectorEntry &entry : element) {
            if (entry.local_coarse_vertex < 0 || entry.local_coarse_vertex >= 3)
                throw std::out_of_range(
                    "Helmholtz corrector local coarse vertex is invalid");
        }
    }

    // setFromTriplets uses the sparse matrix's 32-bit StorageIndex for its
    // prefix sums.  At k=128 the element-wise list contains more than INT_MAX
    // entries because the same fine-row/coarse-column pair is repeated by all
    // incident coarse elements.  The prefix sum then overflows and Eigen writes
    // outside its inner-index array.  Aggregate one coarse column at a time so
    // only unique CSC entries are retained.
    const auto incidence = build_node_incidence(coarse);
    std::vector<ColumnEntries> columns(coarse.nodes.size());
    std::size_t unique_nonzeros = 0;
    constexpr std::size_t storage_index_max =
        static_cast<std::size_t>(
            std::numeric_limits<ComplexSparseMatrix::StorageIndex>::max());
    for (int column = 0; column < static_cast<int>(coarse.nodes.size()); ++column) {
        columns[column] = assemble_column(
            column, fine_node_count, incidence, correctors);
        if (columns[column].size() > storage_index_max - unique_nonzeros) {
            std::ostringstream message;
            message
                << "assembled Helmholtz corrector has more than "
                << storage_index_max
                << " unique nonzeros; rebuild the sparse type with a 64-bit "
                   "StorageIndex";
            throw std::overflow_error(message.str());
        }
        unique_nonzeros += columns[column].size();
    }

    ComplexSparseMatrix matrix(fine_node_count, static_cast<int>(coarse.nodes.size()));
    matrix.reserve(static_cast<Eigen::Index>(unique_nonzeros));
    for (int column = 0; column < static_cast<int>(columns.size()); ++column) {
        matrix.startVec(column);
        for (const ColumnEntry &entry : columns[column])
            matrix.insertBack(entry.row, column) = entry.value;
    }
    matrix.finalize();
    return matrix;
}

ComplexSparseMatrix build_helmholtz_corrected_basis(
    const Eigen::SparseMatrix<double> &coarse_to_fine,
    const TriMesh &coarse,
    int fine_node_count,
    const std::vector<HelmholtzElementCorrector> &correctors) {
    ComplexSparseMatrix basis = coarse_to_fine.cast<Complex>();
    basis -= build_helmholtz_corrector_matrix(
        coarse, fine_node_count, correctors);
    basis.prune(Complex(0.0, 0.0), 1e-14);
    basis.makeCompressed();
    return basis;
}

} // namespace lod2d::helmholtz
