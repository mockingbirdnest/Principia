
#pragma once

#include <type_traits>

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {


// A trait to treat types that have a norm uniformly (using Abs for quantities
// or double, and Norm for multivectors).
template<typename T,
         typename = std::void_t<decltype(quantities::Abs(std::declval<T>()))>>
struct Normed : base::not_constructible {
  using NormType = T;
  static NormType Norm(T const& vector);
};

template<typename T>
struct Normed<T, std::void_t<decltype(std::declval<T>().Norm())>>
    : base::not_constructible {
  using NormType = decltype(std::declval<T>().Norm());
  static NormType Norm(T const& vector);
};

using Instant = Point<quantities::Time>;
template<typename Frame>
using Displacement = Vector<quantities::Length, Frame>;
template<typename Frame>
using Position = Point<Displacement<Frame>>;
template<typename Frame>
using Velocity = Vector<quantities::Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<quantities::AngularFrequency, Frame>;

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, quantities::Length, OrthogonalMap>;

template<typename Frame>
using InertiaTensor =
    SymmetricBilinearForm<quantities::MomentOfInertia, Frame, Bivector>;

}  // namespace geometry
}  // namespace principia

#include "geometry/named_quantities_body.hpp"
