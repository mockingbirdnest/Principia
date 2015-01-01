#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/pair.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Displacement;
using principia::geometry::Pair;
using principia::geometry::Position;
using principia::geometry::Velocity;

namespace principia {
namespace physics {

// This class is analogous to the pair which is its base class, except that it
// exports properly-named selectors.  It is implicitly convertible in both
// directions, so clients can generally ignore the difference.  Note however
// that creating a DegreesOfFreedom involves a copy so clients might want to use
// the base type (probably declared as |auto|) when they don't need to access
// the members.
template<typename Frame>
class DegreesOfFreedom : public Pair<Position<Frame>, Velocity<Frame>> {
 public:
  DegreesOfFreedom(Position<Frame> const& position,
                   Velocity<Frame> const& velocity);

  // Not explicit, the point of this constructor is to convert implicitly.
  DegreesOfFreedom(
      Pair<Position<Frame>, 
           Velocity<Frame>> const& base);  // NOLINT(runtime/explicit)

  Position<Frame> const& position() const;
  Velocity<Frame> const& velocity() const;
};

// This class is analogous to the vector class underlying DegreesOfFreedom,
// except that it exports properly-named selectors.  The same comments as above
// apply.
template<typename Frame>
class RelativeDegreesOfFreedom
    : public Pair<Displacement<Frame>, Velocity<Frame>> {
 public:
  // Not explicit, the point of this constructor is to convert implicitly.
  RelativeDegreesOfFreedom(
      Pair<Displacement<Frame>, 
           Velocity<Frame>> const& base);  // NOLINT(runtime/explicit)

  Displacement<Frame> const& displacement() const;
  Velocity<Frame> const& velocity() const;
};

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights);

}  // namespace physics
}  // namespace principia

#include "physics/degrees_of_freedom_body.hpp"
