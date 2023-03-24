#pragma once

#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/identity.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _homothecy {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::quantities::_named_quantities;

template<typename Scalar, typename FromFrame, typename ToFrame>
class Homothecy : public LinearMap<FromFrame, ToFrame> {
  static_assert(FromFrame::handedness == ToFrame::handedness,
                "Cannot perform an homothecy between frames with different "
                "handedness");

 public:
  Cube<Scalar> Determinant() const;

  Homothecy<Inverse<Scalar>, ToFrame, FromFrame> Inverse() const;

  template<typename VScalar>
  Vector<Product<VScalar, Scalar>, ToFrame> operator()(
      Vector<VScalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename Scalar,
           template<typename, typename> typename Multivector>
  SymmetricBilinearForm<Scalar, ToFrame, Multivector> operator()(
      SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const;

  template<typename T>
  typename base::Mappable<Homothecy, T>::type operator()(T const& t) const;

  template<template<typename, typename> typename LinearMap>
  LinearMap<FromFrame, ToFrame> Forget() const;

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static Homothecy ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Homothecy*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static Homothecy ReadFromMessage(serialization::Homothecy const& message);

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;
};

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
Homothecy<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    Homothecy<LScalar, ThroughFrame, ToFrame> const& left,
    Homothecy<RScalar, FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame, typename Scalar>
std::ostream& operator<<(
    std::ostream& out,
    Homothecy<FromFrame, ToFrame, Scalar> const& homothecy);

}  // namespace internal

using internal::Homothecy;

}  // namespace _homothecy
}  // namespace geometry
}  // namespace principia