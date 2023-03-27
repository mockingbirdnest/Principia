#pragma once

#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/identity.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _conformal_map {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::quantities::_named_quantities;

// This is really a *linear* conformal map, but we don't call it
// LinearConformalMap because that's a mouthful.
template<typename Scalar, typename FromFrame, typename ToFrame>
class ConformalMap : public LinearMap<ConformalMap<Scalar, FromFrame, ToFrame>,
                                      FromFrame, ToFrame> {
 public:
  template<typename S = Scalar,
           typename = std::enable_if_t<!std::is_floating_point_v<S> &&
                                       !std::is_integral_v<S>>>
  ConformalMap(Scalar const& scale,
               OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map);

  Cube<Scalar> Determinant() const;

  ConformalMap<Inverse<Scalar>, ToFrame, FromFrame> Inverse() const;

  template<typename S = Scalar,
           typename = std::enable_if_t<std::is_floating_point_v<S> ||
                                       std::is_integral_v<S>>>
  static ConformalMap Identity();

  template<typename VScalar>
  Vector<Product<VScalar, Scalar>, ToFrame> operator()(
      Vector<VScalar, FromFrame> const& vector) const;

  template<typename T>
  typename base::Mappable<ConformalMap, T>::type operator()(T const& t) const;

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static ConformalMap ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::ConformalMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static ConformalMap ReadFromMessage(
      serialization::ConformalMap const& message);

 private:
  struct PrivateConstructor {};
  ConformalMap(PrivateConstructor,
               Scalar const& scale,
               OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map);

  Scalar const scale_;
  OrthogonalMap<FromFrame, ToFrame> const orthogonal_map_;

  template<typename L, typename R,
           typename From, typename Through, typename To>
  friend ConformalMap<Product<L, R>, From, To> operator*(
      ConformalMap<L, Through, To> const& left,
      ConformalMap<R, From, Through> const& right);

  template<typename S, typename From, typename To>
  friend std::ostream& operator<<(
      std::ostream& out,
      ConformalMap<S, From, To> const& conformal_map);
};

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    ConformalMap<LScalar, ThroughFrame, ToFrame> const& left,
    ConformalMap<RScalar, FromFrame, ThroughFrame> const& right);

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    ConformalMap<Scalar, FromFrame, ToFrame> const& conformal_map);

}  // namespace internal

using internal::ConformalMap;

}  // namespace _conformal_map
}  // namespace geometry
}  // namespace principia

#include "geometry/conformal_map_body.hpp"
