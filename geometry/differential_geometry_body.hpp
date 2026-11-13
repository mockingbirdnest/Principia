#pragma once

#include "geometry/differential_geometry.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace geometry {
namespace _differential_geometry {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_fixed_arrays;

template<typename Scalar, typename S, int s, bool uh>
struct Generator<Scalar, FixedVector<S, s, uh>> {
  using Gradient = FixedVector<Quotient<Scalar, S>, s, uh>;
  static FixedMatrix<double, s, s, uh> InnerProductForm();
};

template<typename Scalar, typename S, typename F>
struct Generator<Scalar, Vector<S, F>> {
  using Gradient = Vector<Quotient<Scalar, S>, F>;
  static SymmetricBilinearForm<double, F, Vector> InnerProductForm();
};

template<typename Scalar, typename V>
struct Generator<Scalar, Point<V>> {
  using Gradient = typename Generator<Scalar, V>::Gradient;
#if _MSC_FULL_VER == 193'632'532 || \
    _MSC_FULL_VER == 193'632'535 || \
    _MSC_FULL_VER == 193'732'822 || \
    _MSC_FULL_VER == 193'833'135
  using InnerProductFormResult =
      decltype(Generator<Scalar, V>::InnerProductForm());
  static InnerProductFormResult InnerProductForm();
#else
  static decltype(Generator<Scalar, V>::InnerProductForm()) InnerProductForm();
#endif
};

template<typename Scalar, typename S, int s, bool uh>
FixedMatrix<double, s, s, uh>
Generator<Scalar, FixedVector<S, s, uh>>::InnerProductForm() {
  FixedMatrix<double, s, s> result{};
  for (int i = 0; i < s; ++i) {
    result(i, i) = 1;
  }
  return result;
}

template<typename Scalar, typename S, typename F>
SymmetricBilinearForm<double, F, Vector>
Generator<Scalar, Vector<S, F>>::InnerProductForm() {
  return geometry::_symmetric_bilinear_form::InnerProductForm<F, Vector>();
}

template<typename Scalar, typename V>
#if _MSC_FULL_VER == 193'632'532 || \
    _MSC_FULL_VER == 193'632'535 || \
    _MSC_FULL_VER == 193'732'822 || \
    _MSC_FULL_VER == 193'833'135
typename Generator<Scalar, Point<V>>::InnerProductFormResult
#else
decltype(Generator<Scalar, V>::InnerProductForm())
#endif
Generator<Scalar, Point<V>>::InnerProductForm() {
  return Generator<Scalar, V>::InnerProductForm();
}

}  // namespace internal
}  // namespace _differential_geometry
}  // namespace geometry
}  // namespace principia
