
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
namespace internal_named_quantities {

using quantities::AngularFrequency;
using quantities::Infinity;
using quantities::Length;
using quantities::MomentOfInertia;
using quantities::Speed;
using quantities::Time;

using Instant = Point<Time>;

CONSTEXPR_INFINITY Instant InfinitePast = Instant() - Infinity<Time>;
CONSTEXPR_INFINITY Instant InfiniteFuture = Instant() + Infinity<Time>;

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

}  // namespace internal_named_quantities

using internal_named_quantities::AngularVelocity;
using internal_named_quantities::Displacement;
using internal_named_quantities::InertiaTensor;
using internal_named_quantities::InfiniteFuture;
using internal_named_quantities::InfinitePast;
using internal_named_quantities::Instant;
using internal_named_quantities::JustAfter;
using internal_named_quantities::JustBefore;
using internal_named_quantities::Position;
using internal_named_quantities::RigidTransformation;
using internal_named_quantities::Velocity;

namespace internal_point {
// We must declare this in the internal namespace where Point is defined so that
// it is found by ADL.
std::string DebugString(const Instant& t);
std::ostream& operator<<(std::ostream& os, const Instant& t);
}  // namespace internal_point

}  // namespace geometry
}  // namespace principia
