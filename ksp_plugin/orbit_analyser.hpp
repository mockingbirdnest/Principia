#pragma once

#include <atomic>
#include <memory>
#include <optional>
#include <thread>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/synchronization/mutex.h"
#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/лидов.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "graphics/colours.hpp"
#include "graphics/graph.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace _orbit_analyser {
namespace internal {

using namespace principia::astronomy::_orbit_ground_track;
using namespace principia::astronomy::_orbit_recurrence;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_лидов;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::graphics::_colours;
using namespace principia::graphics::_graph;
using namespace principia::ksp_plugin::_frames;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_quantities;

// The `OrbitAnalyser` asynchronously integrates a trajectory, and computes
// orbital elements, recurrence, and ground track properties of the resulting
// orbit.
class OrbitAnalyser {
 public:
  class ElementGraphs {
   public:
    struct PlotOptions {
      // The width of all graphs.
      std::int64_t width;
      // The height of time series graphs: a(t), e(t), i(t), ω(t), Ω(t),
      // h_pe(t), h_ap(t).
      std::int64_t time_series_height;
      // The background colour of all graphs.
      RGBA32 background_colour;
      // The axis colour of all graphs.
      RGB24 axis_colour;
      // The colour used for the locus of the eccentricity vector, as well as
      // for the time series of its polar coordinates e(t) and ω(t), and for the
      // lines of constant extremal e in the Лидов graph.
      RGB24 eccentricity_vector_colour;
      // The colour used for the time series of the inclination i(t) and for the
      // lines of constant extremal i in the Лидов graph.
      RGB24 inclination_colour;
      // The colour used for the time series Ω(t).
      RGB24 longitude_of_ascending_node_colour;
      // The colour used for the time series whose dimension is a length: a(t),
      // h_pe(t), h_ap(t).
      RGB24 distance_colour;
      RGB24 лидов_parameter_colour;
      ЛидовGrid лидов_grid;

      friend bool operator==(PlotOptions const& left,
                             PlotOptions const& right) = default;
    };

    // `*elements` must outlive the constructed object.
    ElementGraphs(not_null<OrbitalElements const*> elements,
                  PlotOptions const& options);

    Graph<double, double> const& eccentricity_vector_graph() const;
    Graph<double, double> const& лидов_graph() const;
    Graph<Instant, Length> const& semimajor_axis_graph() const;
    Graph<Instant, double> const& eccentricity_graph() const;
    Graph<Instant, Angle> const& inclination_graph() const;
    Graph<Instant, Angle> const& longitude_of_ascending_node_graph() const;
    Graph<Instant, Angle> const& argument_of_periapsis_graph() const;
    Graph<Instant, Length> const& periapsis_distance_graph() const;
    Graph<Instant, Length> const& apoapsis_distance_graph() const;

    PlotOptions const& plot_options() const;

    void SetЛидовGrid(ЛидовGrid лидов_grid);

   private:
    static not_null<std::unique_ptr<Graph<double, double>>> MakeЛидовGraph(
        OrbitalElements const& elements,
        PlotOptions const& options);

    OrbitalElements const& elements_;
    Graph<double, double> eccentricity_vector_graph_;
    // The Лидов graph can be changed independently of the others by changing
    // the grid options, and Graph is not assignable (the dimensions are fixed
    // at construction), hence the indirection.
    not_null<std::unique_ptr<Graph<double, double>>> лидов_graph_;
    Graph<Instant, Length> semimajor_axis_graph_;
    Graph<Instant, double> eccentricity_graph_;
    Graph<Instant, Angle> inclination_graph_;
    Graph<Instant, Angle> longitude_of_ascending_node_graph_;
    Graph<Instant, Angle> argument_of_periapsis_graph_;
    Graph<Instant, Length> periapsis_distance_graph_;
    Graph<Instant, Length> apoapsis_distance_graph_;
    PlotOptions plot_options_;
  };

  // The analysis stores the computed orbital characteristics.  It is publicly
  // mutable via `SetRecurrence` and `ResetRecurrence` to allow the caller to
  // consider a nominal recurrence other than the one deduced from the orbital
  // elements: analysing the precomputed ground track with respect to a
  // different recurrence is relatively cheap, so it is inconvenient to wait for
  // a whole new analysis to do so, but doing it at every frame is still
  // wasteful, so we cache that in the `Analysis`.
  // It is likewise mutable via `SetPlotOptions`, for the same reasons.
  // Eventually it may make sense to make plotting asynchronous, but it should
  // still likely be separate from analysis.
  class Analysis {
   public:
    Instant const& first_time() const;
    Time const& mission_duration() const;
    RotatingBody<Barycentric> const* primary() const;
    std::optional<Interval<Length>> radial_distance_interval() const;
    std::optional<Instant> first_collision() const;
    std::optional<Instant> first_collision_risk() const;
    std::optional<Instant> first_reentry() const;
    std::optional<OrbitalElements> const& elements() const;
    std::optional<OrbitRecurrence> const& recurrence() const;
    std::optional<OrbitGroundTrack> const& ground_track() const;
    // `equatorial_crossings().has_value()` if and only if
    // `recurrence().has_value && ground_track().has_value()`;
    // `*equatorial_crossings()` is
    //   ground_track()->equator_crossing_longitudes(
    //       *recurrence(), /*first_ascending_pass_index=*/1)
    // precomputed to avoid performing this calculation at every frame.
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes> const&
    equatorial_crossings() const;

    // Sets `recurrence`, updating `equatorial_crossings` if needed.
    void SetRecurrence(OrbitRecurrence const& recurrence);
    // Resets `recurrence` to a value deduced from `*elements` by
    // `OrbitRecurrence::ClosestRecurrence`, or to nullopt if
    // `!elements.has_value()`, updating `equatorial_crossings` if needed.
    void ResetRecurrence();

    // Null if the plot options have not been set.
    ElementGraphs const* element_graphs() const;

    // Recomputes graphs if the options have changed.
    void SetPlotOptions(ElementGraphs::PlotOptions const& options);

   private:
    explicit Analysis(Instant const& first_time);

    Instant first_time_;
    Time mission_duration_;
    RotatingBody<Barycentric> const* primary_ = nullptr;
    std::optional<Interval<Length>> radial_distance_interval_;
    std::optional<Instant> first_collision_;
    std::optional<Instant> first_collision_risk_;
    std::optional<Instant> first_reentry_;
    std::optional<OrbitalElements> elements_;
    std::optional<OrbitRecurrence> closest_recurrence_;
    std::optional<OrbitRecurrence> recurrence_;
    std::optional<OrbitGroundTrack> ground_track_;
    std::optional<OrbitGroundTrack::EquatorCrossingLongitudes>
        equatorial_crossings_;
    std::unique_ptr<ElementGraphs> element_graphs_;

    friend class OrbitAnalyser;
  };

  struct Parameters {
    Instant first_time;
    DegreesOfFreedom<Barycentric> first_degrees_of_freedom;
    Time mission_duration;
    // The analyser may compute the trajectory up to `extended_mission_duration`
    // to ensure that at least one revolution is analysed.
    std::optional<Time> extended_mission_duration;
  };

  OrbitAnalyser(not_null<Ephemeris<Barycentric>*> ephemeris,
                Ephemeris<Barycentric>::FixedStepParameters
                    analysed_trajectory_parameters);

  virtual ~OrbitAnalyser();

  // Cancels any computation in progress, causing the next call to
  // `RequestAnalysis` to be processed as fast as possible.
  void Interrupt();

  // Sets the parameters that will be used for the computation of the next
  // analysis.
  void RequestAnalysis(Parameters const& parameters);

  // The last value passed to `RequestAnalysis`.
  std::optional<Parameters> const& last_parameters() const;

  // Sets `analysis()` to the latest computed analysis.
  void RefreshAnalysis();

  // Mutable so that the caller can call `SetRecurrence` and `ResetRecurrence`.
  Analysis* analysis();

  // The result is in [0, 1]; it tracks the progress of the computation of the
  // next analysis.  Note that a new analysis may be ready even if this is not
  // equal to 1, if the analyser is working on a subsequent request.
  double progress_of_next_analysis() const;

 private:
  using PrimaryCentred = Frame<struct PrimaryCentredTag, NonRotating>;

  // Finds the primary body and analyze our orbit around it.
  absl::Status AnalyseOrbit(Parameters const& parameters);

  // Locates the body with the smallest osculating period and returns it and its
  // period.  This function may be stopped.
  absl::Status FindBodyWithSmallestOsculatingPeriod(
      Parameters const& parameters,
      RotatingBody<Barycentric> const*& primary,
      Time& smallest_osculating_period);

  // Flows the `trajectory` with a fixed step integrator using the given
  // `parameters`.  This is done in small increments and
  // `progress_of_next_analysis_` is updated after each increment to be able to
  // display a progress bar.  This function may be stopped.
  absl::Status FlowWithProgressBar(
      Parameters const& parameters,
      Time const& analysis_duration,
      DiscreteTrajectory<Barycentric>& trajectory);

  // If we can find a sun, computes its mean motion around the primary if it
  // doesn't require too long an integration.  If there is no sun, or the
  // integration would take too long, `mean_sun` is set to `std::nullopt`.
  absl::Status ComputeMeanSunIfPossible(
      Parameters const& parameters,
      BodyCentredNonRotatingReferenceFrame<Barycentric, PrimaryCentred> const&
          primary_centred,
      std::optional<OrbitGroundTrack::MeanSun>& mean_sun);

  // Converts the `trajectory` to the given `primary_centred` frame.  This
  // function may be stopped.
  static absl::StatusOr<DiscreteTrajectory<PrimaryCentred>> ToPrimaryCentred(
      BodyCentredNonRotatingReferenceFrame<Barycentric, PrimaryCentred> const&
          primary_centred,
      DiscreteTrajectory<Barycentric> const& trajectory);

  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  Ephemeris<Barycentric>::FixedStepParameters const
      analysed_trajectory_parameters_;

  std::optional<Parameters> last_parameters_;

  std::optional<Analysis> analysis_;

  mutable absl::Mutex lock_;
  std::jthread analyser_;
  // The `analyser_` is idle:
  // — if it is not joinable, e.g. because it was stopped by `Interrupt()`, or
  // — if it is done computing `next_analysis_` and has stopped or is about to
  //   stop executing.
  // If it is joined once idle (and joinable), it will not attempt to acquire
  // `lock_`.
  bool analyser_idle_ ABSL_GUARDED_BY(lock_) = true;
  // `next_analysis_` is set by the `analyser_` thread; it is read and cleared
  // by the main thread.
  std::optional<Analysis> next_analysis_ ABSL_GUARDED_BY(lock_);
  // `progress_of_next_analysis_` is set by the `analyser_` thread; it tracks
  // progress in computing `next_analysis_`.
  std::atomic<double> progress_of_next_analysis_ = 0;
};

}  // namespace internal

using internal::OrbitAnalyser;

}  // namespace _orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
