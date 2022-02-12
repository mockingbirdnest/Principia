#pragma once

#include "numerics/davenport_q_method.hpp"

#include "geometry/r3x3_matrix.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_davenport_q_method {

using geometry::KroneckerProduct;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::R3x3Matrix;
using quantities::Infinity;

template<typename FromFrame, typename ToFrame, typename Weight>
Rotation<FromFrame, ToFrame> DavenportQMethod(
    std::vector<Vector<double, FromFrame>> const& a,
    std::vector<Vector<double, ToFrame>> const& b,
    std::vector<Weight> const& weights) {
  std::int64_t const size = a.size();
  CHECK_EQ(size, b.size());
  CHECK_EQ(size, weights.size());

  // Compute the attitude profile matrix.
  R3x3Matrix<Weight> B;
  for (int i = 0; i < size; ++i) {
    auto const w_b_ᵗa = weights[i] * KroneckerProduct(b[i].coordinates(),
                                                      a[i].coordinates());
    B += w_b_ᵗa;
  }

  R3Element<Weight> const z{B(1, 2) - B(2, 1),
                            B(2, 0) - B(0, 2),
                            B(0, 1) - B(1, 0)};
  auto const S = B + B.Transpose();
  auto const μ = B.Trace();

  // Form the K matrix.
  FixedMatrix<Weight, 4, 4> K;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      K(i, j) = S(i, j);
      if (i == j) {
        K(i, j) -= μ;
      }
    }
  }
  for (int i = 0; i < 3; ++i) {
    K(3, i) = z[i];
    K(i, 3) = z[i];
  }
  K(3, 3) = μ;

  // Compute its eigensystem.
  auto const eigensystem = ClassicalJacobi(K, /*max_iterations=*/20);
  auto const& eigenvalues = eigensystem.eigenvalues;
  auto const& rotation = eigensystem.rotation;

  // Find the most positive eigenvalue and the corresponding eigenvector;
  auto most_positive_eigenvalue = -Infinity<Weight>;
  Quaternion eigenvector;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    if (eigenvalues[i] > most_positive_eigenvalue) {
      most_positive_eigenvalue = eigenvalues[i];
      eigenvector = Quaternion(
          rotation(3, i),
          R3Element<double>(rotation(0, i), rotation(1, i), rotation(2, i)));
    }
  }

  // The conjugation is because [And17] uses active rotations, but our rotations
  // are passive.
  return Rotation<FromFrame, ToFrame>(eigenvector.Conjugate());
}

}  // namespace internal_davenport_q_method
}  // namespace numerics
}  // namespace principia
