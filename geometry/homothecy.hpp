#pragma once

#include "base/algebra.hpp"
#include "base/concepts.hpp"
#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/space.hpp"
#include "quantities/concepts.hpp"

namespace principia {
namespace geometry {
namespace _homothecy {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_concepts;
using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_space;
using namespace principia::quantities::_concepts;

template<convertible_to_quantity Scalar, typename FromFrame, typename ToFrame>
class Homothecy : public LinearMap<Homothecy<Scalar, FromFrame, ToFrame>,
                                   FromFrame, ToFrame> {
  static_assert(FromFrame::handedness == ToFrame::handedness,
                "Cannot perform an homothecy between frames with different "
                "handedness");

 public:
  explicit Homothecy(Scalar const& scale);

  Cube<Scalar> Determinant() const;

  Homothecy<Inverse<Scalar>, ToFrame, FromFrame> Inverse() const;

  template<quantity VScalar>
  Vector<Product<VScalar, Scalar>, ToFrame> operator()(
      Vector<VScalar, FromFrame> const& vector) const;

  AngularVelocity<ToFrame> operator()(
      AngularVelocity<FromFrame> const& angular_velocity) const;

  template<typename T>
  typename Mappable<Homothecy, T>::type operator()(T const& t) const;

  template<template<typename, typename, typename> typename ConformalMap>
  ConformalMap<Scalar, FromFrame, ToFrame> Forget() const;

  static Homothecy Identity()
    requires(std::is_floating_point_v<Scalar> || std::is_integral_v<Scalar>);

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  static Homothecy ReadFromMessage(serialization::LinearMap const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

  void WriteToMessage(not_null<serialization::Homothecy*> message) const;
  static Homothecy ReadFromMessage(serialization::Homothecy const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

 private:
  Scalar const scale_;

  template<convertible_to_quantity S, typename From, typename To>
  friend class Homothecy;

  template<convertible_to_quantity L, convertible_to_quantity R,
           typename From, typename Through, typename To>
  friend Homothecy<Product<L, R>, From, To> operator*(
      Homothecy<L, Through, To> const& left,
      Homothecy<R, From, Through> const& right);

  template<convertible_to_quantity S, typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Homothecy<S, From, To> const& homothecy);
};

template<convertible_to_quantity LScalar, convertible_to_quantity RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
Homothecy<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    Homothecy<LScalar, ThroughFrame, ToFrame> const& left,
    Homothecy<RScalar, FromFrame, ThroughFrame> const& right);

template<convertible_to_quantity Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    Homothecy<Scalar, FromFrame, ToFrame> const& homothecy);

}  // namespace internal

using internal::Homothecy;

}  // namespace _homothecy
}  // namespace geometry
}  // namespace principia

#include "geometry/homothecy_body.hpp"
