#pragma once

#include "physics/frame.hpp"

#include "geometry/named_quantities.hpp"

using principia::geometry::Position;

namespace principia {
namespace physics {

// Default-initialized to {0, 0, 0}.
template<typename Tag, Tag tag, bool is_inertial>
Position<Frame<Tag, tag, is_inertial>> const
Frame<Tag, tag, is_inertial>::origin;

}  // namespace physics
}  // namespace principia
