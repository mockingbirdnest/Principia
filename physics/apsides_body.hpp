#pragma once

#include "physics/apsides.hpp"

namespace principia {
namespace physics {
namespace internal_apsides {

using geometry::Barycentre;
using geometry::Instant;
using geometry::Sign;
using numerics::Hermite3;
using quantities::Length;
using quantities::Square;
using quantities::Variation;

template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    typename DiscreteTrajectory<Frame>::Iterator const begin,
                    typename DiscreteTrajectory<Frame>::Iterator const end,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides) {
  std::experimental::optional<Instant> previous_time;
  std::experimental::optional<DegreesOfFreedom<Frame>>
      previous_degrees_of_freedom;
  std::experimental::optional<Square<Length>> previous_squared_distance;
  std::experimental::optional<Variation<Square<Length>>>
      previous_squared_distance_derivative;

  for (auto it = begin; it != end; ++it) {
    Instant const time = it.time();
    DegreesOfFreedom<Frame> const degrees_of_freedom = it.degrees_of_freedom();
    DegreesOfFreedom<Frame> const body_degrees_of_freedom =
        reference.EvaluateDegreesOfFreedom(time);
    RelativeDegreesOfFreedom<Frame> const relative =
        degrees_of_freedom - body_degrees_of_freedom;
    Square<Length> const squared_distance =
        InnerProduct(relative.displacement(), relative.displacement());
    // This is the derivative of |squared_distance|.
    Variation<Square<Length>> const squared_distance_derivative =
        2.0 * InnerProduct(relative.displacement(), relative.velocity());

    if (previous_squared_distance_derivative &&
        Sign(squared_distance_derivative) !=
            Sign(*previous_squared_distance_derivative)) {
      CHECK(previous_time &&
            previous_degrees_of_freedom &&
            previous_squared_distance);

      // The derivative of |squared_distance| changed sign.  Construct a Hermite
      // approximation of |squared_distance| and find its extrema.
      Hermite3<Instant, Square<Length>> const
          squared_distance_approximation(
              {*previous_time, time},
              {*previous_squared_distance, squared_distance},
              {*previous_squared_distance_derivative,
               squared_distance_derivative});
      std::set<Instant> const extrema =
          squared_distance_approximation.FindExtrema();

      // Now look at the extrema and check that exactly one is in the required
      // time interval.  This is normally the case, but it can fail due to
      // ill-conditioning.
      Instant apsis_time;
      int valid_extrema = 0;
      for (auto const& extremum : extrema) {
        if (extremum >= *previous_time && extremum <= time) {
          apsis_time = extremum;
          ++valid_extrema;
        }
      }
      if (valid_extrema != 1) {
        // Something went wrong when finding the extrema of
        // |squared_distance_approximation|. Use a linear interpolation of
        // |squared_distance_derivative| instead.
        apsis_time = Barycentre<Instant, Variation<Square<Length>>>(
            {time, *previous_time},
            {*previous_squared_distance_derivative,
             -squared_distance_derivative});
      }

      // Now that we know the time of the apsis, construct a Hermite
      // approximation of the position of the body, and use it to derive its
      // degrees of freedom.  Note that an extremum of
      // |squared_distance_approximation| is in general not an extremum for
      // |position_approximation|: the distance computed using the latter is a
      // 6th-degree polynomial.  However, approximating this polynomial using a
      // 3rd-degree polynomial would yield |squared_distance_approximation|, so
      // we shouldn't be far from the truth.
      Hermite3<Instant, Position<Frame>> position_approximation(
          {*previous_time, time},
          {previous_degrees_of_freedom->position(),
           degrees_of_freedom.position()},
          {previous_degrees_of_freedom->velocity(),
           degrees_of_freedom.velocity()});
      DegreesOfFreedom<Frame> const apsis_degrees_of_freedom(
          position_approximation.Evaluate(apsis_time),
          position_approximation.EvaluateDerivative(apsis_time));
      if (Sign(squared_distance_derivative).Negative()) {
        apoapsides.Append(apsis_time, apsis_degrees_of_freedom);
      } else {
        periapsides.Append(apsis_time, apsis_degrees_of_freedom);
      }
    }

    previous_time = time;
    previous_degrees_of_freedom = degrees_of_freedom;
    previous_squared_distance = squared_distance;
    previous_squared_distance_derivative = squared_distance_derivative;
  }
}

}  // namespace internal_apsides
}  // namespace physics
}  // namespace principia
