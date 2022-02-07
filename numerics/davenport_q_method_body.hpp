#pragma once

#include "numerics/davenport_q_method.hpp"

#include "geometry/r3x3_matrix.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_computations.hpp"

namespace principia {
namespace numerics {
namespace internal_davenport_q_method {

using geometry::KroneckerProduct;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::R3x3Matrix;

template<typename FromFrame, typename ToFrame, typename Weight>
Quaternion DavenportQMethod(std::vector<Vector<double, FromFrame>> const& a,
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

  // Find its eigenvector closest to the identity.
  FixedVector<double, 4> const q₀({0, 0, 0, 1});
  auto const eigensystem = RayleighQuotientIteration(K, q₀);
  auto const& q = eigensystem.eigenvector;

  return Quaternion(q[3], R3Element<double>(q[0], q[1], q[2]));
}

}  // namespace internal_davenport_q_method
}  // namespace numerics
}  // namespace principia
