#pragma once
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory_types {
namespace internal {

template<typename Frame>
value_type<Frame>::value_type(Instant const& time,
                              DegreesOfFreedom<Frame> const& degrees_of_freedom)
    : time(time), degrees_of_freedom(degrees_of_freedom) {}

template<typename Frame>
template<std::size_t i, typename Self>
constexpr auto&& value_type<Frame>::get(this Self&& self) {
  if constexpr (i == 0) {
    return self.time;
  } else if constexpr (i == 1) {
    return self.degrees_of_freedom;
  } else {
    static_assert(i < 2, "Index out of bounds in get<value_type>");
  }
}

template<typename Frame>
bool Earlier::operator()(value_type<Frame> const& left,
                         value_type<Frame> const& right) const {
  return left.time < right.time;
}

template<typename Frame>
bool Earlier::operator()(Instant const& left,
                         value_type<Frame> const& right) const {
  return left < right.time;
}

template<typename Frame>
bool Earlier::operator()(value_type<Frame> const& left,
                         Instant const& right) const {
  return left.time < right;
}

}  // namespace internal
}  // namespace _discrete_trajectory_types
}  // namespace physics
}  // namespace principia
