#pragma once

#include "astronomy/orbit_ground_track.hpp"

#include <limits>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/space.hpp"
#include "numerics/angle_reduction.hpp"
#include "physics/apsides.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace astronomy {
namespace _orbit_ground_track {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_space;
using namespace principia::numerics::_angle_reduction;
using namespace principia::physics::_apsides;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

// Note that the origin of this celestial longitude is arbitrary: it is not the
// node of the orbit around the sun (the equinox).  If `PrimaryCentred` is the
// GCRS, this is the right ascension (with respect to the mean equinox of
// J2000).
template<typename PrimaryCentred>
Angle CelestialLongitude(Position<PrimaryCentred> const& q) {
  // TODO(egg): `ToSpherical` is wasteful, as we discard the latitude.
  return (q - PrimaryCentred::origin).coordinates().ToSpherical().longitude;
}

// The resulting angles are neither normalized nor unwound.
template<typename PrimaryCentred, typename Inertial>
std::vector<Angle> PlanetocentricLongitudes(
    DiscreteTrajectory<PrimaryCentred> const& nodes,
    RotatingBody<Inertial> const& primary) {
  std::vector<Angle> longitudes;
  longitudes.reserve(nodes.size());
  for (auto const& [time, degrees_of_freedom] : nodes) {
    longitudes.push_back(
        CelestialLongitude(degrees_of_freedom.position()) -
        primary.AngleAt(time) - π / 2 * Radian);
  }
  return longitudes;
}

// The resulting angle is not normalized.
template<typename Iterator>
Angle MeanSolarTime(Iterator const& it,
                    OrbitGroundTrack::MeanSun const& mean_sun) {
  Time const t = it->time - mean_sun.epoch;
  return π * Radian + CelestialLongitude(it->degrees_of_freedom.position()) -
         (mean_sun.mean_longitude_at_epoch +
          (2 * π * Radian * t / mean_sun.year));
}

template<typename PrimaryCentred>
Interval<Angle> MeanSolarTimesOfNodes(
    DiscreteTrajectory<PrimaryCentred> const& nodes,
    OrbitGroundTrack::MeanSun const& mean_sun) {
  Interval<Angle> mean_solar_times;
  std::optional<Angle> mean_solar_time;
  for (auto node = nodes.begin(); node != nodes.end(); ++node) {
    if (mean_solar_time.has_value()) {
      mean_solar_time =
          UnwindFrom(*mean_solar_time, MeanSolarTime(node, mean_sun));
    } else {
      mean_solar_time = ReduceAngle<0.0, 2 * π>(MeanSolarTime(node, mean_sun));
    }
    mean_solar_times.Include(*mean_solar_time);
  }
  return mean_solar_times;
}

// Returns the interval spanned by the unwound angles
//   longitudes_of_equatorial_crossings[n] - initial_offset - n * Δλᴇ,
// where Δλᴇ is `nominal_recurrence.equatorial_shift()`, and the first angle is
// in [0, 2π].
inline Interval<Angle> ReducedLongitudesOfEquatorialCrossings(
    std::vector<Angle> const& longitudes_of_equatorial_crossings,
    OrbitRecurrence const& nominal_recurrence,
    Angle const& initial_offset) {
  Interval<Angle> reduced_longitudes;
  std::optional<Angle> reduced_longitude;
  for (int n = 0; n != longitudes_of_equatorial_crossings.size(); ++n) {
    Angle const longitude = longitudes_of_equatorial_crossings[n];
    Angle const offset_longitude =
        longitude - initial_offset - n * nominal_recurrence.equatorial_shift();
    reduced_longitude = reduced_longitude.has_value()
                            ? UnwindFrom(*reduced_longitude, offset_longitude)
                            : ReduceAngle<0.0, 2 * π>(offset_longitude);
    reduced_longitudes.Include(*reduced_longitude);
  }
  return reduced_longitudes;
}

inline Interval<Angle>
OrbitGroundTrack::EquatorCrossingLongitudes::longitudes_reduced_to_pass(
    int const pass_index) const {
  // `shift` applies the number of equatorial shifts corresponding to the pass
  // number; `reduction` puts the midpoint of the shifted interval in [0, 2π].
  Angle shift;
  Interval<Angle> longitudes;
  if (pass_index % 2 == 1) {
    shift = ((pass_index - 1) / 2) * nominal_recurrence_.equatorial_shift();
    longitudes = ascending_longitudes_reduced_to_pass_1_;
  } else {
    shift = ((pass_index - 2) / 2) * nominal_recurrence_.equatorial_shift();
    longitudes = descending_longitudes_reduced_to_pass_2_;
  }
  Angle const reduction =
      -std::floor((longitudes.midpoint() + shift) / (2 * π * Radian)) *
      (2 * π * Radian);
  return {longitudes.min + shift + reduction,
          longitudes.max + shift + reduction};
}

inline OrbitGroundTrack::EquatorCrossingLongitudes::EquatorCrossingLongitudes(
    OrbitRecurrence const& nominal_recurrence)
    : nominal_recurrence_(nominal_recurrence) {}

template<typename PrimaryCentred, typename Inertial>
absl::StatusOr<OrbitGroundTrack> OrbitGroundTrack::ForTrajectory(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    RotatingBody<Inertial> const& primary,
    std::optional<MeanSun> const& mean_sun) {
  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  OrbitGroundTrack ground_track;
  RETURN_IF_ERROR(ComputeNodes(trajectory,
                               trajectory.begin(),
                               trajectory.end(),
                               /*t_max=*/InfiniteFuture,
                               Vector<double, PrimaryCentred>({0, 0, 1}),
                               /*max_points=*/std::numeric_limits<int>::max(),
                               ascending_nodes,
                               descending_nodes));
  if (mean_sun.has_value()) {
    if (!ascending_nodes.empty()) {
      ground_track.mean_solar_times_of_ascending_nodes_ =
          MeanSolarTimesOfNodes(ascending_nodes, *mean_sun);
    }
    if (!descending_nodes.empty()) {
      ground_track.mean_solar_times_of_descending_nodes_ =
          MeanSolarTimesOfNodes(descending_nodes, *mean_sun);
    }
  }
  ground_track.longitudes_of_equator_crossings_of_ascending_passes_ =
      PlanetocentricLongitudes(ascending_nodes, primary);
  ground_track.longitudes_of_equator_crossings_of_descending_passes_ =
      PlanetocentricLongitudes(descending_nodes, primary);
  ground_track.first_descending_pass_before_first_ascending_pass_ =
      !ascending_nodes.empty() && !descending_nodes.empty() &&
      descending_nodes.front().time < ascending_nodes.front().time;
  return ground_track;
}

inline OrbitGroundTrack::EquatorCrossingLongitudes
OrbitGroundTrack::equator_crossing_longitudes(
    OrbitRecurrence const& nominal_recurrence,
    int const first_ascending_pass_index) const {
  EquatorCrossingLongitudes equator_crossings{nominal_recurrence};
  if (longitudes_of_equator_crossings_of_ascending_passes_.empty()) {
    return equator_crossings;
  }
  Angle const δ = nominal_recurrence.grid_interval();
  Angle const λ₀ = longitudes_of_equator_crossings_of_ascending_passes_.front();
  Angle const Δλᴇ = nominal_recurrence.equatorial_shift();
  Angle initial_offset =
      std::floor(λ₀ / δ) * δ + ((first_ascending_pass_index - 1) / 2) * Δλᴇ;
  equator_crossings.ascending_longitudes_reduced_to_pass_1_ =
      ReducedLongitudesOfEquatorialCrossings(
          longitudes_of_equator_crossings_of_ascending_passes_,
          nominal_recurrence,
          initial_offset);
  if (first_descending_pass_before_first_ascending_pass_) {
    // Since the first descending pass is before the first ascending pass, using
    // the same offset in `ReducedLongitudesOfEquatorialCrossings` would compute
    // the longitude reduced to pass 0; adjust by one equatorial shift to get
    // pass 2.
    initial_offset -= Δλᴇ;
  }
  equator_crossings.descending_longitudes_reduced_to_pass_2_ =
      ReducedLongitudesOfEquatorialCrossings(
          longitudes_of_equator_crossings_of_descending_passes_,
          nominal_recurrence,
          initial_offset);
  return equator_crossings;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::mean_solar_times_of_ascending_nodes() const {
  return mean_solar_times_of_ascending_nodes_;
}

inline std::optional<Interval<Angle>> const&
OrbitGroundTrack::mean_solar_times_of_descending_nodes() const {
  return mean_solar_times_of_descending_nodes_;
}

}  // namespace internal
}  // namespace _orbit_ground_track
}  // namespace astronomy
}  // namespace principia
