#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Position;
using principia::geometry::Velocity;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

template<typename Frame>
struct DegreesOfFreedom {
  DegreesOfFreedom(Position<Frame> const& position,
                   Velocity<Frame> const& velocity);
  Position<Frame> position;
  Velocity<Frame> velocity;

  template<typename Weight>
  class BarycentreCalculator {
   public:
    BarycentreCalculator() = default;
    ~BarycentreCalculator() = default;

    void Add(DegreesOfFreedom const& degrees_of_freedom, Weight const& weight);
    DegreesOfFreedom const Get() const;

   private:
    bool empty_ = true;
    decltype(std::declval<Position<Frame>() *
             std::declval<Weight>()) position_weighted_sum_;
    decltype(std::declval<Velocity<Frame>() *
             std::declval<Weight>()) velocity_weighted_sum_;
    Weight weight_;

    // We need a reference position to convert points into vectors.  We pick a
    // default constructed Position<> as it doesn't introduce any inaccuracies
    // in the computations.
    static Position<Frame> const reference_position_;
  };
};

template<typename Frame>
bool operator==(DegreesOfFreedom<Frame> const& left,
                DegreesOfFreedom<Frame> const& right);

template<typename Frame>
bool operator!=(DegreesOfFreedom<Frame> const& left,
                DegreesOfFreedom<Frame> const& right);

template<typename Frame, typename Weight>
DegreesOfFreedom<Frame> Barycentre(
    std::vector<DegreesOfFreedom<Frame>> const& degrees_of_freedom,
    std::vector<Weight> const& weights);

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         DegreesOfFreedom<Frame> const& degrees_of_freedom);

}  // namespace physics
}  // namespace principia

#include "physics/degrees_of_freedom_body.hpp"
