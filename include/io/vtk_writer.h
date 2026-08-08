#pragma once

#include "mesh/types.h"

#include <complex>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string_view>

namespace lod2d::io {

struct VtkDoubleFieldView {
    std::string_view name;
    std::span<const double> values;
};

struct VtkComplexFieldView {
    std::string_view name;
    std::span<const std::complex<double>> values;
};

struct VtkIntFieldView {
    std::string_view name;
    std::span<const int> values;
};

struct VtkUInt64FieldView {
    std::string_view name;
    std::span<const std::uint64_t> values;
};

/// Non-owning field views for one triangular mesh.
///
/// Complex fields are written as two Float64 arrays named `<name>_real` and
/// `<name>_imag`. All spans must remain valid until write_vtu returns.
struct VtuDataView {
    std::span<const VtkDoubleFieldView> point_double;
    std::span<const VtkComplexFieldView> point_complex;
    std::span<const VtkIntFieldView> point_int;
    std::span<const VtkUInt64FieldView> point_uint64;

    std::span<const VtkDoubleFieldView> cell_double;
    std::span<const VtkComplexFieldView> cell_complex;
    std::span<const VtkIntFieldView> cell_int;
    std::span<const VtkUInt64FieldView> cell_uint64;
};

/// Stream an ASCII VTK XML UnstructuredGrid (`.vtu`) without copying mesh or
/// field arrays. The writer creates parent directories when needed.
void write_vtu(
    const std::filesystem::path &path,
    const TriMesh &mesh,
    const VtuDataView &data = {});

/// Write explicit physical boundary edges as VTK line cells with an Int32
/// `boundary_tag` cell field (Dirichlet=1, Robin=2).
void write_boundary_vtu(
    const std::filesystem::path &path,
    const TriMesh &mesh);

} // namespace lod2d::io
