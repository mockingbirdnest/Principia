#pragma once

#include "geometry/differential_geometry.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/concepts.hpp"

namespace principia {
namespace geometry {
namespace _differential_geometry {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_concepts;

template<typename Scalar, typename S, int s, bool uh>
struct Generator<Scalar, FixedVector<S, s, uh>> {
  using Gradient = FixedVector<Quotient<Scalar, S>, s, uh>;
  using InverseHessian = FixedMatrix<Quotient<Square<S>, Scalar>, s, s, uh>;
};

template<convertible_to_quantity Scalar, quantity S, typename F>
struct Generator<Scalar, Vector<S, F>> {
  using Gradient = Vector<Quotient<Scalar, S>, F>;
  using InverseHessian =
      VectorSymmetricBilinearForm<Quotient<Square<S>, Scalar>, F>;
};

template<typename Scalar, quantity Q>
struct Generator<Scalar, Q> {
  using Gradient = Quotient<Scalar, Q>;
  using InverseHessian =
      QuantitySymmetricBilinearForm<Quotient<Square<Q>, Scalar>>;
};

template<typename Scalar, typename V>
struct Generator<Scalar, Point<V>> {
  using Gradient = typename Generator<Scalar, V>::Gradient;
  using InverseHessian = typename Generator<Scalar, V>::InverseHessian;
};

}  // namespace internal
}  // namespace _differential_geometry
}  // namespace geometry
}  // namespace principia
