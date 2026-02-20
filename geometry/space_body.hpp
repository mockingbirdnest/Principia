#pragma once

#include "geometry/space.hpp"

#include "geometry/serialization.hpp"

namespace principia {
namespace geometry {
namespace _direct_sum {
namespace internal {

using namespace principia::geometry::_serialization;

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
    not_null<serialization::Pair*> const message) const {
  PointOrMultivectorSerializer<Position<Frame>, serialization::Pair::Element>::
      WriteToMessage(position_, message->mutable_t1());
  PointOrMultivectorSerializer<Velocity<Frame>, serialization::Pair::Element>::
      WriteToMessage(velocity_, message->mutable_t2());
}

template<typename Frame>
DirectSum<Position<Frame>, Velocity<Frame>>
DirectSum<Position<Frame>, Velocity<Frame>>::ReadFromMessage(
    serialization::Pair const& message) {
  auto const position = PointOrMultivectorSerializer<
      Position<Frame>,
      serialization::Pair::Element>::ReadFromMessage(message.t1());
  auto const velocity = PointOrMultivectorSerializer<
      Velocity<Frame>,
      serialization::Pair::Element>::ReadFromMessage(message.t2());
  return DirectSum(position, velocity);
}

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
    not_null<serialization::Pair*> const message) const {
  PointOrMultivectorSerializer<Displacement<Frame>,
                               serialization::Pair::Element>::
      WriteToMessage(displacement_, message->mutable_t1());
  PointOrMultivectorSerializer<Velocity<Frame>, serialization::Pair::Element>::
      WriteToMessage(velocity_, message->mutable_t2());
}

template<typename Frame>
DirectSum<Displacement<Frame>, Velocity<Frame>>
DirectSum<Displacement<Frame>, Velocity<Frame>>::ReadFromMessage(
    serialization::Pair const& message) {
  auto const displacement = PointOrMultivectorSerializer<
      Displacement<Frame>,
      serialization::Pair::Element>::ReadFromMessage(message.t1());
  auto const velocity = PointOrMultivectorSerializer<
      Velocity<Frame>,
      serialization::Pair::Element>::ReadFromMessage(message.t2());
  return DirectSum(displacement, velocity);
}

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

namespace base {
namespace _mappable {
namespace internal {

template<typename Functor, typename Frame>
typename Mappable<Functor,
                  DirectSum<Displacement<Frame>, Velocity<Frame>>>::type
Mappable<Functor, DirectSum<Displacement<Frame>, Velocity<Frame>>>::Do(
    Functor const& functor,
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& direct_sum) {
  return type(functor(direct_sum.displacement()),
              functor(direct_sum.velocity()));
}

}  // namespace internal
}  // namespace _mappable
}  // namespace base

}  // namespace principia
