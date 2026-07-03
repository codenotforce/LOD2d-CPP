#pragma once

#include "mesh/types.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>
#include <functional>

namespace lod2d::helmholtz {

using Complex = std::complex<double>;
using ComplexVector = Eigen::VectorXcd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexSparseMatrix = Eigen::SparseMatrix<Complex>;
using ComplexTriplet = Eigen::Triplet<Complex>;
using ComplexFunction = std::function<Complex(const Point2 &)>;
using ComplexGradientFunction = std::function<Eigen::Vector2cd(const Point2 &)>;

} // namespace lod2d::helmholtz
