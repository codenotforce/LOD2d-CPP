#include "io/vtk_writer.h"
#include "helmholtz/model.h"

#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>

using namespace lod2d;
using namespace lod2d::io;

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

std::string read_text(const std::filesystem::path &path) {
    std::ifstream input(path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

} // namespace

int main() {
    const std::filesystem::path output =
        std::filesystem::temp_directory_path() / "lod2d_vtk_export_test.vtu";
    try {
        TriMesh mesh;
        mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        mesh.elems = {{0, 1, 3}, {1, 2, 3}};

        const std::array<double, 4> scalar{0.0, 1.0, 2.0, 3.0};
        const std::array<std::complex<double>, 4> complex_values{
            std::complex<double>(1.0, -1.0),
            std::complex<double>(2.0, -2.0),
            std::complex<double>(3.0, -3.0),
            std::complex<double>(4.0, -4.0)};
        const std::array<int, 2> levels{2, 3};
        const std::array<std::uint64_t, 2> ids{10, 11};

        const std::array point_double{
            VtkDoubleFieldView{"temperature", scalar}};
        const std::array point_complex{
            VtkComplexFieldView{"u", complex_values}};
        const std::array cell_int{
            VtkIntFieldView{"level", levels}};
        const std::array cell_uint64{
            VtkUInt64FieldView{"element_id", ids}};

        VtuDataView data;
        data.point_double = point_double;
        data.point_complex = point_complex;
        data.cell_int = cell_int;
        data.cell_uint64 = cell_uint64;
        write_vtu(output, mesh, data);

        const std::string xml = read_text(output);
        require(xml.find("NumberOfPoints=\"4\"") != std::string::npos,
                "VTU point count is missing");
        require(xml.find("NumberOfCells=\"2\"") != std::string::npos,
                "VTU cell count is missing");
        require(xml.find("Name=\"u_real\"") != std::string::npos,
                "VTU complex real field is missing");
        require(xml.find("Name=\"u_imag\"") != std::string::npos,
                "VTU complex imaginary field is missing");
        require(xml.find("Name=\"element_id\"") != std::string::npos,
                "VTU UInt64 cell field is missing");
        require(xml.find("0 1 3 1 2 3") != std::string::npos,
                "VTU connectivity is incorrect");

        bool rejected = false;
        try {
            const std::array<double, 3> short_values{1.0, 2.0, 3.0};
            const std::array invalid{
                VtkDoubleFieldView{"short", short_values}};
            VtuDataView invalid_data;
            invalid_data.point_double = invalid;
            write_vtu(output, mesh, invalid_data);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "VTU writer accepted a mismatched point field");

        const std::filesystem::path boundary_output =
            std::filesystem::temp_directory_path() / "lod2d_test_boundary.vtu";
        write_boundary_vtu(
            boundary_output, lod2d::helmholtz::make_helmholtz_l_shape_mesh());
        const std::string boundary_xml = read_text(boundary_output);
        require(boundary_xml.find("Name=\"boundary_tag\"") != std::string::npos,
                "boundary VTU omitted the boundary_tag field");
        require(boundary_xml.find("NumberOfCells=\"8\"") != std::string::npos,
                "boundary VTU has the wrong edge count");
        std::filesystem::remove(boundary_output);

        std::filesystem::remove(output);
        std::cout << "VTU streaming export schema passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::error_code ignored;
        std::filesystem::remove(output, ignored);
        std::cerr << "test_vtk_export failed: " << error.what() << '\n';
        return 1;
    }
}
