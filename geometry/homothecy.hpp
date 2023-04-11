#pragma once

#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/space.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _homothecy {
namespace internal {

using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;

template<typename Scalar, typename FromFrame, typename ToFrame>
class Homothecy : public LinearMap<Homothecy<Scalar, FromFrame, ToFrame>,
                                   FromFrame, ToFrame> {
  static_assert(FromFrame::handedness == ToFrame::handedness,
                "Cannot perform an homothecy between frames with different "
                "handedness");

 public:
  explicit Homothecy(Scalar const& scale);

  Cube<Scalar> Determinant() const;

  Homothecy<Inverse<Scalar>, ToFrame, FromFrame> Inverse() const;

  template<typename VScalar>
  Vector<Product<VScalar, Scalar>, ToFrame> operator()(
      Vector<VScalar, FromFrame> const& vector) const;

  AngularVelocity<ToFrame> operator()(
      AngularVelocity<FromFrame> const& angular_velocity) const;

  template<typename T>
  typename Mappable<Homothecy, T>::type operator()(T const& t) const;

  template<template<typename, typename, typename> typename ConformalMap>
  ConformalMap<Scalar, FromFrame, ToFrame> Forget() const;

  template<typename S = Scalar,
           typename = std::enable_if_t<std::is_floating_point_v<S> ||
                                       std::is_integral_v<S>>>
  static Homothecy Identity();

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<is_serializable_v<F> &&
                                       is_serializable_v<T>>>
  static Homothecy ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Homothecy*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<is_serializable_v<F> &&
                                       is_serializable_v<T>>>
  static Homothecy ReadFromMessage(serialization::Homothecy const& message);

 private:
  Scalar const scale_;

  template<typename S, typename From, typename To>
  friend class Homothecy;

  template<typename L, typename R,
           typename From, typename Through, typename To>
  friend Homothecy<Product<L, R>, From, To> operator*(
      Homothecy<L, Through, To> const& left,
      Homothecy<R, From, Through> const& right);

  template<typename S, typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Homothecy<S, From, To> const& homothecy);
};

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
Homothecy<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    Homothecy<LScalar, ThroughFrame, ToFrame> const& left,
    Homothecy<RScalar, FromFrame, ThroughFrame> const& right);

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    Homothecy<Scalar, FromFrame, ToFrame> const& homothecy);

}  // namespace internal

using internal::Homothecy;

}  // namespace _homothecy
}  // namespace geometry
}  // namespace principia

#include "geometry/homothecy_body.hpp"
