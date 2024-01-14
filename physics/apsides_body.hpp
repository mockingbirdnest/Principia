#pragma once

#include "physics/apsides.hpp"

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
#include "numerics/matrix_computations.hpp"
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
using namespace principia::numerics::_approximation;
using namespace principia::numerics::_hermite3;
using namespace principia::numerics::_matrix_computations;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;

constexpr double max_error_relative_to_radius = 1e-4;

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

  // Initialize the iterators.  After this block |it| designates the first
  // periapsis and |previous_it| designates the previous apoapsis, if there is
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

  std::vector<Interval<Instant>> intervals;
  for (; it != apsides_times.end(); previous_it = it, ++it) {
    // Here |it| designates a periapsis, and |previous_it| the previous
    // apoapsis, if any.
    Instant const periapsis_time = *it;

    // No collision is possible if the periapsis is above |max_radius|.
    if (squared_distance_from_centre(periapsis_time) < max_radius²) {
      Interval<Instant> interval;
      if (previous_it == apsides_times.end()) {
        // No previous periapsis can only happen the first time through the
        // loop.
        CHECK_EQ(periapsis_time, trajectory.t_min());
        interval.min = periapsis_time;
      } else {
        Instant const apoapsis_time = *previous_it;
        CHECK_LE(apoapsis_time, periapsis_time);

        if (squared_distance_from_centre(apoapsis_time) > max_radius²) {
          // The periapsis is below |max_radius| and the preceding apoapsis is
          // above.  Find the intersection point.
          interval.min = Brent(
              [max_radius²,
               &squared_distance_from_centre](Instant const& time) {
                return squared_distance_from_centre(time) - max_radius²;
              },
              apoapsis_time,
              periapsis_time);
        } else {
          // An apoapsis below |max_radius| can only happen the first time
          // through the loop because the nested loop below skips the other
          // ones.
          CHECK_EQ(apoapsis_time, trajectory.t_min());
          interval.min = apoapsis_time;
        }
      }

      // Loop until we find an apoapsis above |max_radius|, or we reach the end
      // of |apsides_time|.  When entering this loop |it| denotes a periapsis
      // and |previous_it| the preceding apoapsis, if any.
      do {
        Instant const periapis_time = *it;
        previous_it = it;
        ++it;
        if (it == apsides_times.end()) {
          interval.max = periapsis_time;
          break;
        }
        // Here |it| designates an apoapsis, and |previous_it| the previous
        // periapsis.
        Instant const apoapsis_time = *it;
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
        }

        // Fill the end of the interval with the current apoapsis time.  This
        // will yield the right value if we reach the end of |apsides_times| in
        // the while below.
        interval.max = apoapsis_time;
        previous_it = it;
        ++it;
      } while (it != apsides_times.end());

      intervals.push_back(interval);

      if (it == apsides_times.end()) {
        break;
      }

      // The outer loop increment will move us to the periapsis right after the
      // apoapsis that we used to compute |interval.max| using Brent.
    } else {
      // Go to the next pair periapsis, apoapsis (in this order).
      previous_it = it;
      ++it;
      if (it == apsides_times.end()) {
        // Reached the end of |apsides_times| with a periapsis above
        // |max_radius|.
        break;
      }
    }
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
    std::function<Length(Angle const& latitude, Angle const& longitude)> const&
        radius) {
  // The frame of the surface of the celestial.
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;

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

  // Interpolate the height above the terrain using a Чебышёв polynomial.
  auto const чебышёв_interpolant = ЧебышёвPolynomialInterpolant(
      height_above_terrain_at_time,
      interval.min,
      interval.max,
      reference_body.mean_radius() * max_error_relative_to_radius);

  auto const companion_matrix = чебышёв_interpolant.FrobeniusCompanionMatrix();
  auto const companion_matrix_schur_decomposition =
      RealSchurDecomposition(companion_matrix, max_error_relative_to_radius);

  auto const& T = companion_matrix_schur_decomposition.T;
  absl::btree_set<double> real_roots;
  for (int i = 0; i < T.rows() - 1; ++i) {
    if (i == T.rows() - 1 && T(i, i - 1) == 0) {
      real_roots.insert(T(i, i));
    } else if (T(i + 1, i) == 0) {
      real_roots.insert(T(i, i));
    }
  }

  for (double r : real_roots) {
    LOG(ERROR)<<r;
  }

  return std::nullopt;
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
