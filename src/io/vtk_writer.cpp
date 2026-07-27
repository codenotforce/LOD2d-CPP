#include "io/vtk_writer.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace lod2d::io {

namespace {

std::string xml_escape(std::string_view text) {
    std::string result;
    result.reserve(text.size());
    for (char ch : text) {
        switch (ch) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '"': result += "&quot;"; break;
        case '\'': result += "&apos;"; break;
        default: result += ch; break;
        }
    }
    return result;
}

void require_field_size(
    std::string_view name,
    std::size_t actual,
    std::size_t expected,
    std::string_view association) {
    if (name.empty()) throw std::invalid_argument("VTU field name must not be empty");
    if (actual != expected) {
        throw std::invalid_argument(
            "VTU " + std::string(association) + " field " + std::string(name)
            + " has " + std::to_string(actual) + " values; expected "
            + std::to_string(expected));
    }
}

template <typename Field>
void validate_fields(
    std::span<const Field> fields,
    std::size_t expected,
    std::string_view association,
    std::unordered_set<std::string> &names) {
    for (const Field &field : fields) {
        require_field_size(field.name, field.values.size(), expected, association);
        if (!names.insert(std::string(field.name)).second) {
            throw std::invalid_argument(
                "duplicate VTU " + std::string(association)
                + " field name: " + std::string(field.name));
        }
    }
}

void validate_complex_names(
    std::span<const VtkComplexFieldView> fields,
    std::unordered_set<std::string> &names,
    std::string_view association) {
    for (const VtkComplexFieldView &field : fields) {
        const std::string real_name = std::string(field.name) + "_real";
        const std::string imag_name = std::string(field.name) + "_imag";
        if (!names.insert(real_name).second || !names.insert(imag_name).second) {
            throw std::invalid_argument(
                "duplicate expanded VTU " + std::string(association)
                + " complex field name: " + std::string(field.name));
        }
    }
}

template <typename Value>
void write_scalar_array(
    std::ostream &output,
    std::string_view vtk_type,
    std::string_view name,
    std::span<const Value> values) {
    output << "        <DataArray type=\"" << vtk_type << "\" Name=\""
           << xml_escape(name) << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (const Value &value : values) output << value << ' ';
    output << "\n        </DataArray>\n";
}

void write_complex_arrays(
    std::ostream &output,
    const VtkComplexFieldView &field) {
    output << "        <DataArray type=\"Float64\" Name=\""
           << xml_escape(std::string(field.name) + "_real")
           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (const std::complex<double> &value : field.values) output << value.real() << ' ';
    output << "\n        </DataArray>\n";

    output << "        <DataArray type=\"Float64\" Name=\""
           << xml_escape(std::string(field.name) + "_imag")
           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    for (const std::complex<double> &value : field.values) output << value.imag() << ' ';
    output << "\n        </DataArray>\n";
}

void write_fields(
    std::ostream &output,
    std::span<const VtkDoubleFieldView> doubles,
    std::span<const VtkComplexFieldView> complexes,
    std::span<const VtkIntFieldView> integers,
    std::span<const VtkUInt64FieldView> uint64_values) {
    for (const VtkDoubleFieldView &field : doubles)
        write_scalar_array(output, "Float64", field.name, field.values);
    for (const VtkComplexFieldView &field : complexes)
        write_complex_arrays(output, field);
    for (const VtkIntFieldView &field : integers)
        write_scalar_array(output, "Int32", field.name, field.values);
    for (const VtkUInt64FieldView &field : uint64_values)
        write_scalar_array(output, "UInt64", field.name, field.values);
}

void validate_mesh(const TriMesh &mesh) {
    const int node_count = static_cast<int>(mesh.nodes.size());
    for (const Triangle &triangle : mesh.elems) {
        for (int node : triangle) {
            if (node < 0 || node >= node_count)
                throw std::invalid_argument("VTU triangle contains an invalid node index");
        }
    }
}

} // namespace

void write_vtu(
    const std::filesystem::path &path,
    const TriMesh &mesh,
    const VtuDataView &data) {
    validate_mesh(mesh);
    const std::size_t point_count = mesh.nodes.size();
    const std::size_t cell_count = mesh.elems.size();

    std::unordered_set<std::string> point_names;
    validate_fields(data.point_double, point_count, "point", point_names);
    validate_fields(data.point_int, point_count, "point", point_names);
    validate_fields(data.point_uint64, point_count, "point", point_names);
    for (const VtkComplexFieldView &field : data.point_complex)
        require_field_size(field.name, field.values.size(), point_count, "point");
    validate_complex_names(data.point_complex, point_names, "point");

    std::unordered_set<std::string> cell_names;
    validate_fields(data.cell_double, cell_count, "cell", cell_names);
    validate_fields(data.cell_int, cell_count, "cell", cell_names);
    validate_fields(data.cell_uint64, cell_count, "cell", cell_names);
    for (const VtkComplexFieldView &field : data.cell_complex)
        require_field_size(field.name, field.values.size(), cell_count, "cell");
    validate_complex_names(data.cell_complex, cell_names, "cell");

    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open VTU output: " + path.string());
    output << std::setprecision(std::numeric_limits<double>::max_digits10);

    output << "<?xml version=\"1.0\"?>\n"
              "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n"
              "  <UnstructuredGrid>\n"
              "    <Piece NumberOfPoints=\"" << point_count
           << "\" NumberOfCells=\"" << cell_count << "\">\n"
              "      <PointData>\n";
    write_fields(
        output,
        data.point_double,
        data.point_complex,
        data.point_int,
        data.point_uint64);
    output << "      </PointData>\n"
              "      <CellData>\n";
    write_fields(
        output,
        data.cell_double,
        data.cell_complex,
        data.cell_int,
        data.cell_uint64);
    output << "      </CellData>\n"
              "      <Points>\n"
              "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const Point2 &point : mesh.nodes)
        output << point.x() << ' ' << point.y() << " 0 ";
    output << "\n        </DataArray>\n"
              "      </Points>\n"
              "      <Cells>\n"
              "        <DataArray type=\"Int32\" Name=\"connectivity\" "
              "format=\"ascii\">\n          ";
    for (const Triangle &triangle : mesh.elems)
        output << triangle[0] << ' ' << triangle[1] << ' ' << triangle[2] << ' ';
    output << "\n        </DataArray>\n"
              "        <DataArray type=\"Int32\" Name=\"offsets\" "
              "format=\"ascii\">\n          ";
    for (std::size_t cell = 0; cell < cell_count; ++cell)
        output << 3 * (cell + 1) << ' ';
    output << "\n        </DataArray>\n"
              "        <DataArray type=\"UInt8\" Name=\"types\" "
              "format=\"ascii\">\n          ";
    for (std::size_t cell = 0; cell < cell_count; ++cell) output << "5 ";
    output << "\n        </DataArray>\n"
              "      </Cells>\n"
              "    </Piece>\n"
              "  </UnstructuredGrid>\n"
              "</VTKFile>\n";
    if (!output) throw std::runtime_error("failed while writing VTU output: " + path.string());
}

} // namespace lod2d::io
