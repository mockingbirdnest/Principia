#pragma once

#include "physics/frame.hpp"

#include "geometry/named_quantities.hpp"

using principia::geometry::Position;

namespace principia {
namespace physics {

// Default-initialized to {0, 0, 0}.
template<int tag, bool is_inertial>
Position<Frame<tag, is_inertial>> const Frame<tag, is_inertial>::origin;

template<int tag, bool is_inertial>
bool const Frame<tag, is_inertial>::is_inertial;

}  // namespace physics
}  // namespace principia
