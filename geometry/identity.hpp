#pragma once

#include "base/concepts.hpp"
#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

FORWARD_DECLARE(
    TEMPLATE(typename Scalar,
             typename Frame,
             template<typename, typename> typename Multivector) class,
    SymmetricBilinearForm,
    FROM(symmetric_bilinear_form),
    INTO(identity));

namespace _identity {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_sign;

// The identity map.
template<typename FromFrame, typename ToFrame>
class Identity : public LinearMap<Identity<FromFrame, ToFrame>,
                                  FromFrame, ToFrame> {
  static_assert(FromFrame::handedness == ToFrame::handedness,
                "Cannot identity frames with different handedness");

 public:
  Identity() = default;

  Sign Determinant() const;

  Identity<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

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
  typename Mappable<Identity, T>::type operator()(T const& t) const;

  template<template<typename, typename> typename LinearMap>
  LinearMap<FromFrame, ToFrame> Forget() const;
  template<template<typename, typename, typename> typename ConformalMap>
  ConformalMap<double, FromFrame, ToFrame> Forget() const;

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  static Identity ReadFromMessage(serialization::LinearMap const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

  void WriteToMessage(not_null<serialization::Identity*> message) const;
  static Identity ReadFromMessage(serialization::Identity const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Identity<FromFrame, ToFrame> operator*(
    Identity<ThroughFrame, ToFrame> const& left,
    Identity<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Identity<FromFrame, ToFrame> const& identity);

}  // namespace internal

using internal::Identity;

}  // namespace _identity
}  // namespace geometry
}  // namespace principia

#include "geometry/identity_body.hpp"
