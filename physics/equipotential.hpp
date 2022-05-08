#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/ephemeris.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using geometry::Bivector;
using geometry::Instant;
using geometry::Position;
using physics::Ephemeris;

//TODO(phl): Avoid Instant?

template<typename Frame>
std::vector<Position<Frame>> ComputeEquipotential(
    Ephemeris<Frame> const& ephemeris,
    Bivector<double, Frame> const& plane,
    Position<Frame> const& position,
    Instant const& t);

}  // namespace internal_equipotential

using internal_equipotential::ComputeEquipotential;

}  // namespace physics
}  // namespace principia

#include "physics/equipotential_body.hpp"
