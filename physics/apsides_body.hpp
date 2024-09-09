#pragma once

#include "physics/apsides.hpp"

#include <algorithm>
#include <list>
#include <optional>
#include <vector>

#include "base/array.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "geometry/space.hpp"
#include "numerics/approximation.hpp"
#include "numerics/hermite3.hpp"
#include "numerics/root_finders.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _apsides {
namespace internal {

using namespace principia::base::_array;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_space;
using namespace principia::numerics::_approximation;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

// This number was selected in game to obtain good performance without missing
// collisions.  It was found that that one evaluation of
// `height_above_terrain_at_time` costs around 15 µs, and that the number of
// evaluations grows roughly linearly with the degree.  This is probably because
// the terrain is very bumpy, so even a high degree doesn't approximate it well,
// and the number of subdivisions end up being largely independent from the
// degree.  It would be faster to use degree 4, but since the height would only
// be evaluated 5 times per interval we could easily miss significant features
// of the terrain.
constexpr int max_чебышёв_degree = 8;
// The error bound `max_collision_error` is guaranteed to be met if the vessel
// is slower than this.
constexpr Speed max_collision_speed = 10'000 * Metre / Second;

template<typename Frame>
void ComputeApsides(Trajectory<Frame> const& reference,
                    Trajectory<Frame> const& trajectory,
                    typename DiscreteTrajectory<Frame>::iterator const begin,
                    typename DiscreteTrajectory<Frame>::iterator const end,
                    Instant const& t_max,
                    int const max_points,
                    DiscreteTrajectory<Frame>& apoapsides,
                    DiscreteTrajectory<Frame>& periapsides) {
  std::optional<Instant> previous_time;
  std::optional<DegreesOfFreedom<Frame>> previous_degrees_of_freedom;
  std::optional<Square<Length>> previous_squared_distance;
  std::optional<Variation<Square<Length>>>
      previous_squared_distance_derivative;

  Instant const effective_t_min = reference.t_min();
  Instant const effective_t_max = std::min(t_max, reference.t_max());
  for (auto it = begin; it != end; ++it) {
    auto const& [time, degrees_of_freedom] = *it;
    if (time < effective_t_min) {
      continue;
    }
    if (time > effective_t_max) {
      break;
    }
    DegreesOfFreedom<Frame> const body_degrees_of_freedom =
        reference.EvaluateDegreesOfFreedom(time);
    RelativeDegreesOfFreedom<Frame> const relative =
        degrees_of_freedom - body_degrees_of_freedom;
    Square<Length> const squared_distance = relative.displacement().Norm²();
    // This is the derivative of `squared_distance`.
    Variation<Square<Length>> const squared_distance_derivative =
        2.0 * InnerProduct(relative.displacement(), relative.velocity());

    if (previous_squared_distance_derivative &&
        Sign(squared_distance_derivative) !=
            Sign(*previous_squared_distance_derivative)) {
      CHECK(previous_time &&
            previous_degrees_of_freedom &&
            previous_squared_distance);

      // The derivative of `squared_distance` changed sign.  Construct a Hermite
      // approximation of `squared_distance` and find its extrema.
      Hermite3<Square<Length>, Instant> const
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
        // `squared_distance_approximation`. Use a linear interpolation of
        // `squared_distance_derivative` instead.
        apsis_time = Barycentre({time, *previous_time},
                                {*previous_squared_distance_derivative,
                                 -squared_distance_derivative});
        // This happens if `*previous_squared_distance_derivative` is 0.  It is
        // indicative of ill-conditioning, but not in and of itself a problem.
        LOG_IF(WARNING, apsis_time == previous_time)
            << "Suspicious apsis at the beginning of a time interval: "
            << apsis_time;
      }

      // This can happen for instance if the square distance is stationary.
      // Safer to give up.
      if (!IsFinite(apsis_time - Instant{})) {
        break;
      }

      // Now that we know the time of the apsis, use a Hermite approximation to
      // derive its degrees of freedom.  Note that an extremum of
      // `squared_distance_approximation` is in general not an extremum for
      // `position_approximation`: the distance computed using the latter is a
      // 6th-degree polynomial.  However, approximating this polynomial using a
      // 3rd-degree polynomial would yield `squared_distance_approximation`, so
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

  // Construct the set of all extrema times (apsides and extremities).  In this
  // set, apoapsides and periapsides alternate.
  absl::btree_set<Instant> apsides_times{trajectory.t_min(),
                                         trajectory.t_max()};
  for (auto const& [time, _] : apoapsides) {
    apsides_times.insert(time);
  }
  for (auto const& [time, _] : periapsides) {
    apsides_times.insert(time);
  }
  if (apsides_times.size() < 2) {
    return {};
  }

  // Initialize the iterators.  After this block `it` designates the first
  // periapsis and `previous_it` designates the previous apoapsis, if there is
  // one, or is past the end.
  absl::btree_set<Instant>::const_iterator it;
  absl::btree_set<Instant>::const_iterator previous_it;
  {
    auto const first_it = apsides_times.begin();
    auto const second_it = std::next(first_it);
    if (squared_distance_from_centre(*first_it) <
        squared_distance_from_centre(*second_it)) {
      previous_it = apsides_times.end();
      it = first_it;
    } else {
      previous_it = first_it;
      it = second_it;
    }
  }

  // Verify that the apoapsides and periapsides alternate in distance from the
  // centre.  This is a critical assumption made by the algorithm below.  If
  // they don't, just give up.
  {
    bool at_apoapsis = false;
    bool at_periapsis = false;
    bool is_anomalous = false;

    auto it = apsides_times.begin();
    Square<Length> distance = squared_distance_from_centre(*it);
    Square<Length> previous_distance;

    for (++it; it != apsides_times.end(); ++it) {
      previous_distance = distance;
      distance = squared_distance_from_centre(*it);
      if (previous_distance < distance) {
        if (at_apoapsis) {
          is_anomalous = true;
          break;
        }
        at_apoapsis = true;
        at_periapsis = false;
      } else if (distance < previous_distance) {
        if (at_periapsis) {
          is_anomalous = true;
          break;
        }
        at_apoapsis = false;
        at_periapsis = true;
      } else {  // distance == previous_distance
        is_anomalous = true;
        break;
      }
    }
    if (is_anomalous) {
      LOG(WARNING) << "Anomalous apsides, not computing collisions";
      return {};
    }
  }

  std::vector<Interval<Instant>> intervals;
  for (; it != apsides_times.end(); previous_it = it, ++it) {
    // Here `it` designates a periapsis, and `previous_it` the previous
    // apoapsis, or is past-the-end if there is no previous apoapsis (only the
    // first time through the outer loop).
    Instant const periapsis_time = *it;

    // No collision is possible if the periapsis is above `max_radius`.
    if (squared_distance_from_centre(periapsis_time) < max_radius²) {
      Interval<Instant> interval;
      if (previous_it == apsides_times.end()) {
        // No previous apoapsis can only happen the first time through the outer
        // loop.
        CHECK_EQ(periapsis_time, trajectory.t_min());
        interval.min = periapsis_time;
      } else {
        Instant const apoapsis_time = *previous_it;
        CHECK_LE(apoapsis_time, periapsis_time);

        if (squared_distance_from_centre(apoapsis_time) > max_radius²) {
          // The periapsis is below `max_radius` and the preceding apoapsis is
          // above.  Find the intersection point.
          interval.min = Brent(
              [max_radius²,
               &squared_distance_from_centre](Instant const& time) {
                return squared_distance_from_centre(time) - max_radius²;
              },
              apoapsis_time,
              periapsis_time);
        } else {
          // An apoapsis below `max_radius` can only happen the first time
          // through the loop because the nested loop below skips the other
          // ones.
          CHECK_EQ(apoapsis_time, trajectory.t_min());
          interval.min = apoapsis_time;
        }
      }
      // Here `interval.min` has been computed and is the time where the
      // trajectory descends below `max_radius`.

      // Loop until we find an apoapsis above `max_radius`, or we reach the end
      // of `apsides_time`.
      do {
        // Here `it` denotes a periapsis below `max_radius` and `previous_it`
        // the preceding apoapsis or is past-the-end if there is no previous
        // apoapsis (only the first time through the outer and inner loops).
        Instant const periapsis_time = *it;
        previous_it = it;
        ++it;

        // Here `it` denotes an apoapsis or is past-the-end.  `previous_it` is
        // the previous periapsis, which is below `max_radius`.
        if (it == apsides_times.end()) {
          break;
        }
        Instant const apoapsis_time = *it;
        CHECK_LE(periapsis_time, apoapsis_time);

        if (squared_distance_from_centre(apoapsis_time) > max_radius²) {
          // The periapsis is below `max_radius` and the following apoapsis is
          // above.  Find the intersection point.
          interval.max = Brent(
              [max_radius²,
               &squared_distance_from_centre](Instant const& time) {
                return squared_distance_from_centre(time) - max_radius²;
              },
              periapsis_time,
              apoapsis_time);
          break;
        }

        previous_it = it;
        ++it;

        // Here `it` denotes a periapsis or is past-the end.  `previous_it` is
        // the previous apoapsis, which is below `max_radius`.  Therefore, the
        // periapsis denoted by `it` is below `max_radius`.
      } while (it != apsides_times.end());

      if (it == apsides_times.end()) {
        // If we reach the end of the set, the last element, denoted by
        // `previous_it`, may be a periapsis (first break in the loop above) or
        // an apoapsis (while condition).  In both cases, `previous_it` denotes
        // the extremity of the trajectory.
        interval.max = *previous_it;
        intervals.push_back(interval);
        break;
      }

      // `interval.max` has been computed using Brent (second break in the loop
      // above).
      intervals.push_back(interval);

      // Here `it` designates an apoapsis above `max_radius`.  `previous_it` is
      // the previous periapsis, which is below `max_radius`.  The outer loop
      // increment moves `previous_it` to denote the apoapsis and `it` to denote
      // the next periapsis or to be past-the-end.
    } else {
      // Here `it` designates a periapsis above `max_radius`, and `previous_it`
      // the previous apoapsis, or is past-the-end if there is no previous
      // apoapsis (only the first time through the outer loop).
      previous_it = it;
      ++it;
      // Here `it` denotes an apoapsis or is past-the-end.  `previous_it` is
      // the previous periapsis, which is above `max_radius`.  Therefore, the
      // apoapsis denoted by `it` is above `max_radius`.
      if (it == apsides_times.end()) {
        // Reached the end of `apsides_times` with a periapsis above
        // `max_radius`.  No interval to produce.
        break;
      }
    }
  }

  // Check the consistency of the result.
  for (std::int64_t i = 1; i < intervals.size(); ++i) {
    CHECK_LE(intervals[i - 1].max, intervals[i].min);
    CHECK_LE(intervals[i].min, intervals[i].max);
  }

  return intervals;
}

template<typename Frame>
std::optional<typename DiscreteTrajectory<Frame>::value_type>
ComputeFirstCollision(
    RotatingBody<Frame> const& reference_body,
    Trajectory<Frame> const& reference,
    Trajectory<Frame> const& trajectory,
    Interval<Instant> const& interval,
    Length const& max_collision_error,
    std::function<Length(Angle const& latitude, Angle const& longitude)> const&
        radius) {
  // The frame of the surface of the celestial.
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;

  std::int64_t number_of_evaluations = 0;

  auto height_above_terrain_at_time = [&number_of_evaluations,
                                       &radius,
                                       &reference,
                                       &reference_body,
                                       &trajectory](Instant const& t) {
    ++number_of_evaluations;
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

  // Subdivide the interpolant if it could have real roots given the current
  // error estimate.
  SubdivisionPredicate<Length, Instant> const subdivide =
      [](auto const& interpolant, Length const& error_estimate) -> bool {
    return interpolant.MayHaveRealRoots(error_estimate);
  };

  // Stop if the interpolant has a real root.  This is the first collision.
  std::optional<Instant> first_collision_time;
  TerminationPredicate<Length, Instant> const stop =
      [&first_collision_time,
       max_collision_error,
       &number_of_evaluations](auto interpolant) -> bool {
    if (interpolant->MayHaveRealRoots()) {
      // The relative error on the roots is choosen so that it corresponds to an
      // absolute error in distance similar to `max_collision_error`, assuming
      // that the speed is below `max_collision_speed`.  We don't care too much
      // about the performance of this computation, because the zero-free test
      // is very efficient.
      auto const& real_roots = interpolant->RealRoots(
          max_collision_error /
          ((interpolant->upper_bound() - interpolant->lower_bound()) *
           max_collision_speed));
      if (real_roots.empty()) {
        VLOG(1) << "No real roots over [" << interpolant->lower_bound() << ", "
                << interpolant->upper_bound() << "]";
        return false;
      } else {
        // The smallest root is the first collision.
        first_collision_time = *real_roots.begin();
        VLOG(1) << "First collision time is " << first_collision_time.value()
                << " with " << number_of_evaluations << " evaluations";
        return true;
      }
    } else {
      VLOG(1) << "No roots over [" << interpolant->lower_bound() << ", "
              << interpolant->upper_bound() << "]";
      return false;
    }
  };

  // Interpolate the height above the terrain using Чебышёв polynomials.
  StreamingAdaptiveЧебышёвPolynomialInterpolant<max_чебышёв_degree>(
      height_above_terrain_at_time,
      interval.min,
      interval.max,
      max_collision_error,
      subdivide,
      stop);
  if (first_collision_time.has_value()) {
    // The first collision.
    return typename DiscreteTrajectory<Frame>::value_type(
        *first_collision_time,
        trajectory.EvaluateDegreesOfFreedom(*first_collision_time));
  } else {
    // No collision.
    return std::nullopt;
  }
}

template<typename Frame, typename Predicate>
absl::Status ComputeNodes(
    Trajectory<Frame> const& trajectory,
    typename DiscreteTrajectory<Frame>::iterator const begin,
    typename DiscreteTrajectory<Frame>::iterator const end,
    Instant const& t_max,
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
    if (time > t_max) {
      break;
    }
    Length const z =
        (degrees_of_freedom.position() - Frame::origin).coordinates().z;
    Speed const z_speed = degrees_of_freedom.velocity().coordinates().z;

    if (previous_z && Sign(z) != Sign(*previous_z)) {
      CHECK(previous_time && previous_z_speed);

      // `z` changed sign.  Construct a Hermite approximation of `z` and find
      // its zeros.
      Hermite3<Length, Instant> const z_approximation(
          {*previous_time, time},
          {*previous_z, z},
          {*previous_z_speed, z_speed});

      Instant node_time;
      if (Sign(z_approximation.Evaluate(*previous_time)) ==
          Sign(z_approximation.Evaluate(time))) {
        // The Hermite approximation is poorly conditioned, let's use a linear
        // approximation
        node_time = Barycentre({*previous_time, time}, {z, -*previous_z});
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
          // `north` is up and we are going up, or `north` is down and we are
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
