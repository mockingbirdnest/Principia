#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/pair.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Displacement;
using principia::geometry::Position;
using principia::geometry::Velocity;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

template<typename Frame>
class DegreesOfFreedom
    : public principia::geometry::Pair<Position<Frame>, Velocity<Frame>> {
 public:
  DegreesOfFreedom(Position<Frame> const& position,
                   Velocity<Frame> const& velocity);

  Position<Frame> const& position() const;
  Velocity<Frame> const& velocity() const;

  //TODO(phl): Move to Pair.
  template<typename Weight>
  class BarycentreCalculator {
   public:
    BarycentreCalculator() = default;
    ~BarycentreCalculator() = default;

    void Add(DegreesOfFreedom const& degrees_of_freedom, Weight const& weight);
    DegreesOfFreedom const Get() const;

   private:
    bool empty_ = true;
    decltype(std::declval<Displacement<Frame>>() *
             std::declval<Weight>()) displacements_weighted_sum_;
    decltype(std::declval<Velocity<Frame>>() *
             std::declval<Weight>()) velocities_weighted_sum_;
    Weight weight_;

    // We need a reference position to convert points into vectors.  We pick a
    // default constructed Position<> as it doesn't introduce any inaccuracies
    // in the computations.
    static Position<Frame> const reference_position_;
  };
};

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights);

}  // namespace physics
}  // namespace principia

#include "physics/degrees_of_freedom_body.hpp"
