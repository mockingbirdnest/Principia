#pragma once

#include <vector>

#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
DegreesOfFreedom<Frame>::DegreesOfFreedom(Position<Frame> const& position,
                                          Velocity<Frame> const& velocity)
    : Pair<Position<Frame>, Velocity<Frame>>(position, velocity) {}

template<typename Frame>
DegreesOfFreedom<Frame>::DegreesOfFreedom(
    Pair<Position<Frame>, Velocity<Frame>> const& base)
    : Pair<Position<Frame>, Velocity<Frame>>(base) {}

template<typename Frame>
DegreesOfFreedom<Frame> DegreesOfFreedom<Frame>::ReadFromMessage(
    serialization::Pair const& message) {
  return Pair<Position<Frame>, Velocity<Frame>>::ReadFromMessage(message);
}

template<typename Frame>
Position<Frame> const& DegreesOfFreedom<Frame>::position() const {
  return this->t1_;
}

template<typename Frame>
Velocity<Frame> const& DegreesOfFreedom<Frame>::velocity() const {
  return this->t2_;
}

template<typename Frame>
RelativeDegreesOfFreedom<Frame>::RelativeDegreesOfFreedom(
    Displacement<Frame> const& displacement,
    Velocity<Frame> const& velocity)
    : Pair<Displacement<Frame>, Velocity<Frame>>(displacement, velocity) {}

template<typename Frame>
RelativeDegreesOfFreedom<Frame>::RelativeDegreesOfFreedom(
    Pair<Displacement<Frame>, Velocity<Frame>> const& base)
    : Pair<Displacement<Frame>, Velocity<Frame>>(base) {}

template<typename Frame>
Displacement<Frame> const&
RelativeDegreesOfFreedom<Frame>::displacement() const {
  return this->t1_;
}

template<typename Frame>
Velocity<Frame> const& RelativeDegreesOfFreedom<Frame>::velocity() const {
  return this->t2_;
}

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights) {
  CHECK_EQ(degrees_of_freedom.size(), weights.size())
      << "Degrees of freedom and weights of unequal sizes";
  CHECK(!degrees_of_freedom.empty()) << "Empty input";
  typename DegreesOfFreedom<Frame>::
      template BarycentreCalculator<Weight> calculator;
  for (size_t i = 0; i < degrees_of_freedom.size(); ++i) {
    calculator.Add(degrees_of_freedom[i], weights[i]);
  }
  return calculator.Get();
  CHECK_EQ(degrees_of_freedom.size(), weights.size());
  CHECK(!degrees_of_freedom.empty());
}

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  out << "{" << degrees_of_freedom.position() << ", "
      << degrees_of_freedom.velocity() << "}";
  return out;
}

template<typename Frame>
std::ostream& operator<<(
    std::ostream& out,
    RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom) {
  out << "{" << relative_degrees_of_freedom.displacement() << ", "
      << relative_degrees_of_freedom.velocity() << "}";
  return out;
}

}  // namespace physics

namespace base {

template<typename Functor, typename Frame>
typename Mappable<Functor, physics::RelativeDegreesOfFreedom<Frame>>::type
Mappable<Functor, physics::RelativeDegreesOfFreedom<Frame>>::Do(
    Functor const& functor,
    physics::RelativeDegreesOfFreedom<Frame> const& relative) {
  return type(functor(relative.displacement()), functor(relative.velocity()));
}

}  // namespace base
}  // namespace principia
