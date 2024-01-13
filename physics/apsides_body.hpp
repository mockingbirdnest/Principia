#pragma once

#include "physics/apsides.hpp"

#include <optional>

#include "base/array.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "geometry/space.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/root_finders.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _apsides {
namespace internal {

using namespace principia::base::_array;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_space;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    Trajectory<Frame> const& trajectory,
                    typename DiscreteTrajectory<Frame>::iterator const begin,
                    typename DiscreteTrajectory<Frame>::iterator const end,
                    int const max_points,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides) {
  std::optional<Instant> previous_time;
  std::optional<DegreesOfFreedom<Frame>> previous_degrees_of_freedom;
  std::optional<Square<Length>> previous_squared_distance;
  std::optional<Variation<Square<Length>>>
      previous_squared_distance_derivative;

  Instant const t_min = reference.t_min();
  Instant const t_max = reference.t_max();
  for (auto it = begin; it != end; ++it) {
    auto const& [time, degrees_of_freedom] = *it;
    if (time < t_min) {
      continue;
    }
    if (time > t_max) {
      break;
    }
    DegreesOfFreedom<Frame> const body_degrees_of_freedom =
        reference.EvaluateDegreesOfFreedom(time);
    RelativeDegreesOfFreedom<Frame> const relative =
        degrees_of_freedom - body_degrees_of_freedom;
    Square<Length> const squared_distance = relative.displacement().Norm²();
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
      BoundedArray<Instant, 2> const extrema =
          squared_distance_approximation.FindExtrema();

      // Now look at the extrema and check that exactly one is in the required
      // time interval.  This is normally the case, but it can fail due to
      // ill-conditioning.
      Instant apsis_time;
      int valid_extrema = 0;
      for (auto const& extremum : extrema) {
        if (extremum > *previous_time && extremum <= time) {
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

      // This can happen for instance if the square distance is stationary.
      // Safer to give up.
      if (!IsFinite(apsis_time - Instant{})) {
        break;
      }

      // Now that we know the time of the apsis, use a Hermite approximation to
      // derive its degrees of freedom.  Note that an extremum of
      // |squared_distance_approximation| is in general not an extremum for
      // |position_approximation|: the distance computed using the latter is a
      // 6th-degree polynomial.  However, approximating this polynomial using a
      // 3rd-degree polynomial would yield |squared_distance_approximation|, so
      // we shouldn't be far from the truth.
      DegreesOfFreedom<Frame> const apsis_degrees_of_freedom =
          trajectory.EvaluateDegreesOfFreedom(apsis_time);
      if (Sign(squared_distance_derivative).is_negative()) {
        apoapsides.Append(apsis_time, apsis_degrees_of_freedom).IgnoreError();
      } else {
        periapsides.Append(apsis_time, apsis_degrees_of_freedom).IgnoreError();
      }
      if (apoapsides.size() >= max_points && periapsides.size() >= max_points) {
        break;
      }
    }

    previous_time = time;
    previous_degrees_of_freedom = degrees_of_freedom;
    previous_squared_distance = squared_distance;
    previous_squared_distance_derivative = squared_distance_derivative;
  }
}

template<typename Frame>
std::vector<Interval<Instant>> ComputeCollisionIntervals(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    DiscreteTrajectory<Frame> const& apoapsides,
    DiscreteTrajectory<Frame> const& periapsides) {
  auto squared_distance_from_centre = [&reference,
                                       &trajectory](Instant const& time) {
    Position<Frame> const body_position =
        reference.EvaluatePosition(time);
    Position<Frame> const trajectory_position =
        trajectory.EvaluatePosition(time);
    Displacement<Frame> const displacement =
        trajectory_position - body_position;
    return displacement.Norm²();
  };

  auto const max_radius² = Pow<2>(reference_body.max_radius());

  // NOTE(phl): This function doesn't consider |t_min| and |t_max| as potential
  // periapsis, only as potential apoapsis.  Unsure if that's right.

  // Construct the list of times of the apoapsides.  This list includes the
  // times |t_min| and |t_max| if they are "on the way" to an apoapsis.
  std::list<Instant> apoapsides_times;
  for (auto const& [time, _] : apoapsides) {
    apoapsides_times.push_back(time);
  }
  if (apoapsides_times.empty() ||
      apoapsides_times.front() > periapsides.begin()->time) {
    apoapsides_times.push_front(trajectory.t_min());
  }
  if (periapsides.empty() ||
      apoapsides_times.back() < std::prev(periapsides.end())->time) {
    apoapsides_times.push_back(trajectory.t_max());
  }
  CHECK_EQ(apoapsides_times.size(), periapsides.size() + 1);

  std::vector<Interval<Instant>> intervals;

  auto ait = apoapsides_times.begin();
  for (auto pit = periapsides.begin(); pit != periapsides.end();) {
    Instant apoapsis_time = *ait;
    Instant periapsis_time = pit->time;
    CHECK_LE(apoapsis_time, periapsis_time);

    // No collision is possible if the periapsis is above |max_radius|.
    if (squared_distance_from_centre(periapsis_time) < max_radius²) {
      Interval<Instant> interval;
      if (squared_distance_from_centre(apoapsis_time) > max_radius²) {
        // The periapsis is below |max_radius| and the preceding apoapsis is
        // above.  Find the intersection point.
        interval.min = Brent(
            [max_radius², &squared_distance_from_centre](Instant const& time) {
              return squared_distance_from_centre(time) - max_radius²;
            },
            apoapsis_time,
            periapsis_time);
      } else {
        // A periapsis below |max_radius| can only happen at the beginning of
        // the list because the loop below skips the other ones.
        CHECK(pit == periapsides.begin());
        interval.min = periapsis_time;
      }

      // Loop until we find an apoapsis above |max_radius|, or we reach the end
      // of |apoapsides_time|.  We only move |periapsis_time| to limit the
      // range that the Brent algorithm has to search.
      for (++ait; ait != apoapsides_times.end(); ++ait, ++pit) {
        CHECK(pit != periapsides.end());
        apoapsis_time = *ait;
        periapsis_time = pit->time;
        CHECK_LE(periapsis_time, apoapsis_time);

        if (squared_distance_from_centre(apoapsis_time) > max_radius²) {
          // The periapsis is below |max_radius| and the following apoapsis is
          // above.  Find the intersection point.
          interval.max = Brent(
              [max_radius²,
               &squared_distance_from_centre](Instant const& time) {
                return squared_distance_from_centre(time) - max_radius²;
              },
              periapsis_time,
              apoapsis_time);
          break;
        } else {
          // Fill the end of the interval with the current apoapsis time.  This
          // will yield the right value if we reach the end of
          // |apoapsides_time|.
          interval.max = apoapsis_time;
        }
      }

      intervals.push_back(interval);

      // Go to the periapsis right after the apoapsis that we used to compute
      // |interval.max|.
      ++pit;
    } else {
      // Go to the next pair periapsis, apoapsis (in this order).
      ++ait;
      ++pit;
    }
  }

  return intervals;
}


template<typename Frame>
typename DiscreteTrajectory<Frame>::value_type ComputeCollision(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    Instant const& first_time,
    Instant const& last_time,
    std::function<Length(Angle const& latitude,
                         Angle const& longitude)> const& radius) {
  // The frame of the surface of the celestial.
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;

  auto squared_distance_from_centre = [&reference,
                                       &trajectory](Instant const& time) {
    Position<Frame> const body_position =
        reference.EvaluatePosition(time);
    Position<Frame> const trajectory_position =
        trajectory.EvaluatePosition(time);
    Displacement<Frame> const displacement =
        trajectory_position - body_position;
    return displacement.Norm²();
  };

  auto const max_radius² = Pow<2>(reference_body.max_radius());
  auto const min_radius² = Pow<2>(reference_body.min_radius());

  // Find the time at which the |trajectory| crosses |max_radius|, if any.
  // We'll start the search from there.  This is cheap because it doesn't call
  // C#.
  Instant max_radius_time;
  if (squared_distance_from_centre(first_time) > max_radius²) {
    max_radius_time = Brent(
        [max_radius², &squared_distance_from_centre](Instant const& time) {
          return squared_distance_from_centre(time) - max_radius²;
        },
        first_time,
        last_time);
  } else {
    max_radius_time = first_time;
  }

  // Similarly, find the time at which the |trajectory| crosses |min_radius|, if
  // any.
  Instant min_radius_time;
  if (squared_distance_from_centre(last_time) < min_radius²) {
    min_radius_time = Brent(
        [min_radius², &squared_distance_from_centre](Instant const& time) {
          return squared_distance_from_centre(time) - min_radius²;
        },
        max_radius_time,
        last_time);
  } else {
    min_radius_time = last_time;
  }

  auto height_above_terrain_at_time = [&radius,
                                       &reference,
                                       &reference_body,
                                       &trajectory](Instant const& t) {
    auto const reference_position = reference.EvaluatePosition(t);
    auto const trajectory_position = trajectory.EvaluatePosition(t);
    Displacement<Frame> const displacement_in_frame =
        trajectory_position - reference_position;

    auto const to_surface_frame =
        reference_body.template ToSurfaceFrame<SurfaceFrame>(t);
    Displacement<SurfaceFrame> const displacement_in_surface =
        to_surface_frame(displacement_in_frame);

    SphericalCoordinates<Length> const spherical_coordinates =
        displacement_in_surface.coordinates().ToSpherical();

    return spherical_coordinates.radius -
           radius(spherical_coordinates.latitude,
                  spherical_coordinates.longitude);
  };

  CHECK_LE(Length{}, height_above_terrain_at_time(max_radius_time));
  CHECK_LE(height_above_terrain_at_time(min_radius_time), Length{});

  Instant const collision_time = Brent(height_above_terrain_at_time,
                                       max_radius_time,
                                       min_radius_time);

  return typename DiscreteTrajectory<Frame>::value_type(
      collision_time,
      trajectory.EvaluateDegreesOfFreedom(collision_time));
}

template<typename Frame, typename Predicate>
absl::Status ComputeNodes(
    Trajectory<Frame> const& trajectory,
    typename DiscreteTrajectory<Frame>::iterator const begin,
    typename DiscreteTrajectory<Frame>::iterator const end,
    Vector<double, Frame> const& north,
    int const max_points,
    DiscreteTrajectory<Frame>& ascending,
    DiscreteTrajectory<Frame>& descending,
    Predicate predicate) {
  static_assert(
      std::is_convertible<decltype(predicate(
                              std::declval<DegreesOfFreedom<Frame>>())),
                          bool>::value,
      "|predicate| must be a predicate on |DegreesOfFreedom<Frame>|");

  std::optional<Instant> previous_time;
  std::optional<Length> previous_z;
  std::optional<Speed> previous_z_speed;

  for (auto it = begin; it != end; ++it) {
    RETURN_IF_STOPPED;
    auto const& [time, degrees_of_freedom] = *it;
    Length const z =
        (degrees_of_freedom.position() - Frame::origin).coordinates().z;
    Speed const z_speed = degrees_of_freedom.velocity().coordinates().z;

    if (previous_z && Sign(z) != Sign(*previous_z)) {
      CHECK(previous_time && previous_z_speed);

      // |z| changed sign.  Construct a Hermite approximation of |z| and find
      // its zeros.
      Hermite3<Instant, Length> const z_approximation(
          {*previous_time, time},
          {*previous_z, z},
          {*previous_z_speed, z_speed});

      Instant node_time;
      if (Sign(z_approximation.Evaluate(*previous_time)) ==
          Sign(z_approximation.Evaluate(time))) {
        // The Hermite approximation is poorly conditioned, let's use a linear
        // approximation
        node_time = Barycentre<Instant, Length>({*previous_time, time},
                                                {z, -*previous_z});
      } else {
        // The normal case, find the intersection with z = 0 using Brent's
        // method.
        node_time = Brent(
            [&z_approximation](Instant const& t) {
              return z_approximation.Evaluate(t);
            },
            *previous_time,
            time);
      }

      DegreesOfFreedom<Frame> const node_degrees_of_freedom =
          trajectory.EvaluateDegreesOfFreedom(node_time);
      if (predicate(node_degrees_of_freedom)) {
        if (Sign(InnerProduct(north, Vector<double, Frame>({0, 0, 1}))) ==
            Sign(z_speed)) {
          // |north| is up and we are going up, or |north| is down and we are
          // going down.
          ascending.Append(node_time, node_degrees_of_freedom).IgnoreError();
        } else {
          descending.Append(node_time, node_degrees_of_freedom).IgnoreError();
        }
        if (ascending.size() >= max_points && descending.size() >= max_points) {
          break;
        }
      }
    }

    previous_time = time;
    previous_z = z;
    previous_z_speed = z_speed;
  }
  return absl::OkStatus();
}

}  // namespace internal
}  // namespace _apsides
}  // namespace physics
}  // namespace principia
