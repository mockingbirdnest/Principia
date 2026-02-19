#pragma once

#include "geometry/space.hpp"

namespace principia {
namespace geometry {
namespace _direct_sum {
namespace internal {

template<typename Frame>
constexpr DirectSum<Position<Frame>, Velocity<Frame>>::DirectSum(
    uninitialized_t) {}

template<typename Frame>
constexpr DirectSum<Position<Frame>, Velocity<Frame>>::DirectSum(
    Position<Frame> const& position,
    Velocity<Frame> const& velocity)
    : position_(position), velocity_(velocity) {}

template<typename Frame>
DirectSum<Position<Frame>, Velocity<Frame>>&
DirectSum<Position<Frame>, Velocity<Frame>>::operator+=(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& right) {
  *this = *this + right;
  return *this;
}

template<typename Frame>
DirectSum<Position<Frame>, Velocity<Frame>>&
DirectSum<Position<Frame>, Velocity<Frame>>::operator-=(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& right) {
  *this = *this - right;
  return *this;
}

template<typename Frame>
template<typename Self>
constexpr auto&& DirectSum<Position<Frame>, Velocity<Frame>>::position(
    this Self&& self) {
  return self.position_;
}
template<typename Frame>
template<typename Self>
constexpr auto&& DirectSum<Position<Frame>, Velocity<Frame>>::velocity(
    this Self&& self) {
  return self.velocity_;
}

template<typename Frame>
void DirectSum<Position<Frame>, Velocity<Frame>>::WriteToMessage(
    not_null<serialization::Pair*> const message) const {}

template<typename Frame>
DirectSum<Position<Frame>, Velocity<Frame>>
DirectSum<Position<Frame>, Velocity<Frame>>::ReadFromMessage(
    serialization::Pair const& message) {}

template<typename Frame>
constexpr DirectSum<Displacement<Frame>, Velocity<Frame>>::DirectSum(
    uninitialized_t) {}

template<typename Frame>
constexpr DirectSum<Displacement<Frame>, Velocity<Frame>>::DirectSum(
    Displacement<Frame> const& displacement,
    Velocity<Frame> const& velocity)
    : displacement_(displacement), velocity_(velocity) {}

template<typename Frame>
DirectSum<Displacement<Frame>, Velocity<Frame>>&
DirectSum<Displacement<Frame>, Velocity<Frame>>::operator+=(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& right) {
  *this = *this + right;
  return *this;
}

template<typename Frame>
DirectSum<Displacement<Frame>, Velocity<Frame>>&
DirectSum<Displacement<Frame>, Velocity<Frame>>::operator-=(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& right) {
  *this = *this - right;
  return *this;
}

template<typename Frame>
template<typename Self>
constexpr auto&& DirectSum<Displacement<Frame>, Velocity<Frame>>::displacement(
    this Self&& self) {
  return self.displacement_;
}
template<typename Frame>
template<typename Self>
constexpr auto&& DirectSum<Displacement<Frame>, Velocity<Frame>>::velocity(
    this Self&& self) {
  return self.velocity_;
}

template<typename Frame>
void DirectSum<Displacement<Frame>, Velocity<Frame>>::WriteToMessage(
    not_null<serialization::Pair*> const message) const {}

template<typename Frame>
DirectSum<Displacement<Frame>, Velocity<Frame>>
DirectSum<Displacement<Frame>, Velocity<Frame>>::ReadFromMessage(
    serialization::Pair const& message) {}

template<std::size_t i, typename Frame>
constexpr auto const& get(
    DirectSum<Position<Frame>, Velocity<Frame>> const& direct_sum) {
  if constexpr (i == 0) {
    return direct_sum.position();
  } else if constexpr (i == 1) {
    return direct_sum.velocity();
  } else {
    static_assert(i < 2,
                  "Index out of bounds in get<DirectSum<Position, Velocity>>");
  }
}

template<std::size_t i, typename Frame>
constexpr auto const& get(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& direct_sum) {
  if constexpr (i == 0) {
    return direct_sum.displacement();
  } else if constexpr (i == 1) {
    return direct_sum.velocity();
  } else {
    static_assert(
        i < 2, "Index out of bounds in get<DirectSum<Displacement, Velocity>>");
  }
}

}  // namespace internal
}  // namespace _direct_sum
}  // namespace geometry
}  // namespace principia
