#pragma once

#include <limits>
#include <string>
#include <type_traits>

#include "base/macros.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

FORWARD_DECLARE_FROM_NEW(orthogonal_map,
                         TEMPLATE(typename FromFrame, typename ToFrame) class,
                         OrthogonalMap);

namespace _named_quantities {
namespace internal {

using quantities::AngularFrequency;
using quantities::Infinity;
using quantities::Length;
using quantities::MomentOfInertia;
using quantities::Speed;
using quantities::Time;
using namespace principia::geometry::_orthogonal_map;

using Instant = Point<Time>;

constexpr Instant InfinitePast = Instant() - Infinity<Time>;
constexpr Instant InfiniteFuture = Instant() + Infinity<Time>;

template<typename Frame>
using Displacement = Vector<Length, Frame>;

template<typename Frame>
using Position = Point<Displacement<Frame>>;

template<typename Frame>
using Velocity = Vector<Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<AngularFrequency, Frame>;

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, Length, OrthogonalMap>;

template<typename Frame>
using InertiaTensor = SymmetricBilinearForm<MomentOfInertia, Frame, Bivector>;

// IEEE 754:2008 nextUp and nextDown for Instants.
// We would like to avoid the terms “up” and “down” when referring to the
// passage of time.  We avoid the term “next” in one direction because of the
// confusability with |std::nextafter|, which has different semantics, and in
// the other because of the awkwardness of the phrase “next before”.
// Defined inline for want of a way to alias functions in C++.
constexpr Instant JustAfter(Instant const t) { return NextUp(t); }
constexpr Instant JustBefore(Instant const t) { return NextDown(t); }

}  // namespace internal

using internal::AngularVelocity;
using internal::Displacement;
using internal::InertiaTensor;
using internal::InfiniteFuture;
using internal::InfinitePast;
using internal::Instant;
using internal::JustAfter;
using internal::JustBefore;
using internal::Position;
using internal::RigidTransformation;
using internal::Velocity;

namespace _named_quantities {
namespace internal {
// We must declare this in the internal namespace where Point is defined so that
// it is found by ADL.
std::string DebugString(const Instant& t);
std::ostream& operator<<(std::ostream& os, const Instant& t);
}  // namespace internal

}  // namespace _named_quantities
}  // namespace _named_quantities
}  // namespace geometry
}  // namespace principia

namespace principia::geometry {
using namespace principia::geometry::_named_quantities;
}  // namespace principia::geometry
