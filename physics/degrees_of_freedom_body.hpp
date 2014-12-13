#pragma once

#include <vector>

#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
DegreesOfFreedom<Frame>::DegreesOfFreedom(Position<Frame> const& position,
                                          Velocity<Frame> const& velocity)
    : position(position),
      velocity(velocity) {}

template<typename Frame>
template<typename Weight>
void DegreesOfFreedom<Frame>::BarycentreCalculator<Weight>::Add(
    DegreesOfFreedom const& degrees_of_freedom,
    Weight const& weight) {
  auto const displacements_weighted_sum_diff =
      (degrees_of_freedom.position - reference_position_) * weight;
  auto const velocity_weighted_sum_diff =
      degrees_of_freedom.velocity * weight;
  if (empty_) {
    displacements_weighted_sum_ = displacements_weighted_sum_diff;
    velocities_weighted_sum_ = velocity_weighted_sum_diff;
    weight_ = weight;
    empty_ = false;
  } else {
    displacements_weighted_sum_ += displacements_weighted_sum_diff;
    velocities_weighted_sum_ += velocity_weighted_sum_diff;
    weight_ += weight;
  }
}

template<typename Frame>
template<typename Weight>
DegreesOfFreedom<Frame> const
DegreesOfFreedom<Frame>::BarycentreCalculator<Weight>::Get() const {
  CHECK(!empty_) << "Empty BarycentreCalculator";
  return {reference_position_ +
              Displacement<Frame>(displacements_weighted_sum_ / weight_),
          Velocity<Frame>(velocities_weighted_sum_ / weight_)};
}

template<typename Frame>
bool operator==(DegreesOfFreedom<Frame> const& left,
                DegreesOfFreedom<Frame> const& right) {
  return left.position == right.position &&
         left.velocity == right.velocity;
}

template<typename Frame>
bool operator!=(DegreesOfFreedom<Frame> const& left,
                DegreesOfFreedom<Frame> const& right) {
  return left.position != right.position ||
         left.velocity != right.velocity;
}

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights) {
  CHECK_EQ(degrees_of_freedom.size(), weights.size())
      << "Degrees of freedom and weights of unequal sizes";
  CHECK(!degrees_of_freedom.empty()) << "Empty input";
  DegreesOfFreedom<Frame>::BarycentreCalculator<Weight> calculator;
  for (size_t i = 0; i < degrees_of_freedom.size(); ++i) {
    calculator.Add(degrees_of_freedom[i], weights[i]);
  }
  return calculator.Get();
  CHECK_EQ(degrees_of_freedom.size(), weights.size());
  CHECK(!degrees_of_freedom.empty());
  // We need a reference position to convert points into vectors.  We pick a
  // default constructed Position<> as it doesn't introduce any inaccuracies in
  // the computations below.
  Position<Frame> const reference_position;
  auto positions_weighted_sum =
      (degrees_of_freedom[0].position -
       reference_position).coordinates() * weights[0];
  auto velocities_weighted_sum =
      degrees_of_freedom[0].velocity.coordinates() * weights[0];
  Weight weight = weights[0];
  for (size_t i = 1; i < degrees_of_freedom.size(); ++i) {
    positions_weighted_sum +=
        (degrees_of_freedom[i].position -
         reference_position).coordinates() * weights[i];
    velocities_weighted_sum +=
        degrees_of_freedom[i].velocity.coordinates() * weights[i];
    weight += weights[i];
  }
  return {reference_position +
              geometry::Displacement<Frame>(positions_weighted_sum / weight),
          geometry::Velocity<Frame>(velocities_weighted_sum / weight)};
}

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  return out << "{" << degrees_of_freedom.position << ", "
                    << degrees_of_freedom.velocity << "}";
}

template<typename Frame>
template<typename Weight>
Position<Frame> const
DegreesOfFreedom<Frame>::BarycentreCalculator<Weight>::reference_position_;

}  // namespace physics
}  // namespace principia
