#pragma once

#include <memory>
#include <string>
#include <typeindex>
#include <type_traits>
#include <utility>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/push_pull_callback.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rp2_point.hpp"
#include "geometry/space.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/orbit_analyser.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/planetarium.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/renderer.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace interface {

// This is used for interfacing, and should only appear in C++ code in tests
// and generated code; we allow ourselves to pollute the `interface` namespace
// with convenience `using`s.

using namespace principia::base::_not_null;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;
using namespace principia::base::_push_pull_callback;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_rp2_point;
using namespace principia::geometry::_space;
using namespace principia::ksp_plugin::_flight_plan;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_iterators;
using namespace principia::ksp_plugin::_orbit_analyser;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_renderer;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Takes ownership of `**pointer` and returns it to the caller.  Nulls
// `*pointer`.  `pointer` must not be null.  No transfer of ownership of
// `*pointer`.
template<typename T>
std::unique_ptr<T> TakeOwnership(T** pointer);
template<typename T>
std::unique_ptr<T[]> TakeOwnershipArray(T** pointer);

#include "ksp_plugin/interface.generated.h"

extern "C" PRINCIPIA_DLL
void __cdecl principia__ActivatePlayer();

extern "C" PRINCIPIA_DLL
void __cdecl principia__ActivateRecorder(bool activate);

extern "C" PRINCIPIA_DLL
void __cdecl principia__InitGoogleLogging();

bool operator==(AdaptiveStepParameters const& left,
                AdaptiveStepParameters const& right);
bool operator==(Burn const& left, Burn const& right);
bool operator==(EquatorialCrossings const& left,
                EquatorialCrossings const& right);
bool operator==(FlightPlanAdaptiveStepParameters const& left,
                FlightPlanAdaptiveStepParameters const& right);
bool operator==(Interval const& left, Interval const& right);
bool operator==(NavigationFrameParameters const& left,
                NavigationFrameParameters const& right);
bool operator==(NavigationManoeuvre const& left,
                NavigationManoeuvre const& right);
bool operator==(NavigationManoeuvreFrenetTrihedron const& left,
                NavigationManoeuvreFrenetTrihedron const& right);
bool operator==(Node const& left, Node const& right);
bool operator==(OrbitAnalysis const& left, OrbitAnalysis const& right);
bool operator==(OrbitRecurrence const& left, OrbitRecurrence const& right);
bool operator==(OrbitalElements const& left, OrbitalElements const& right);
bool operator==(QP const& left, QP const& right);
bool operator==(QPRW const& left, QPRW const& right);
bool operator==(SolarTimesOfNodes const& left, SolarTimesOfNodes const& right);
bool operator==(TQP const& left, TQP const& right);
bool operator==(WXYZ const& left, WXYZ const& right);
bool operator==(XY const& left, XY const& right);
bool operator==(XYZ const& left, XYZ const& right);

// Conversions between interchange data and typed data.

physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters
FromAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters);
std::pair<physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters,
          physics::_ephemeris::Ephemeris<
              Barycentric>::GeneralizedAdaptiveStepParameters>
FromFlightPlanAdaptiveStepParameters(FlightPlanAdaptiveStepParameters const&
                                         flight_plan_adaptive_step_parameters);

Renderer::Node FromNode(Plugin const& plugin,
                        Node const& node);

template<typename T>
T FromQP(QP const& qp);
template<>
DegreesOfFreedom<World> FromQP<DegreesOfFreedom<World>>(QP const& qp);
template<>
RelativeDegreesOfFreedom<AliceSun>
FromQP<RelativeDegreesOfFreedom<AliceSun>>(QP const& qp);
template<>
RelativeDegreesOfFreedom<World>
FromQP<RelativeDegreesOfFreedom<World>>(QP const& qp);

Quaternion FromWXYZ(WXYZ const& wxyz);

R3Element<double> FromXYZ(XYZ const& xyz);
template<typename T>
T FromXYZ(XYZ const& xyz);
template<>
Position<World> FromXYZ<Position<World>>(XYZ const& xyz);
template<>
Position<EccentricPart> FromXYZ<Position<EccentricPart>>(XYZ const& xyz);
template<>
Velocity<Frenet<NavigationFrame>>
FromXYZ<Velocity<Frenet<NavigationFrame>>>(XYZ const& xyz);
template<>
AngularVelocity<World> FromXYZ(XYZ const& xyz);
template<>
R3Element<MomentOfInertia> FromXYZ<R3Element<MomentOfInertia>>(XYZ const& xyz);

AdaptiveStepParameters ToAdaptiveStepParameters(
    physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters);
FlightPlanAdaptiveStepParameters ToFlightPlanAdaptiveStepParameters(
    physics::_ephemeris::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    physics::_ephemeris::Ephemeris<
        Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters);

KeplerianElements ToKeplerianElements(
    physics::_kepler_orbit::KeplerianElements<Barycentric> const&
        keplerian_elements);

Node ToNode(Plugin const& plugin,
            Renderer::Node const& node);

QP ToQP(DegreesOfFreedom<World> const& dof);
QP ToQP(RelativeDegreesOfFreedom<AliceSun> const& relative_dof);

// Ownership of the status and its message is transferred to the caller.
Status* ToNewStatus(absl::Status const& status);

WXYZ ToWXYZ(Quaternion const& quaternion);

XY ToXY(RP2Point<Length, Camera> const& rp2_point);

XYZ ToXYZ(R3Element<double> const& r3_element);
XYZ ToXYZ(Position<World> const& position);
XYZ ToXYZ(Vector<double, World> const& direction);
XYZ ToXYZ(Velocity<Frenet<NavigationFrame>> const& velocity);
XYZ ToXYZ(Bivector<AngularMomentum, World> const& angular_momentum);

template<typename T>
Interval ToInterval(geometry::_interval::Interval<T> const& interval);

// Conversions between interchange data and typed data that depend on the state
// of the plugin.
Instant FromGameTime(Plugin const& plugin, double t);
double ToGameTime(Plugin const& plugin, Instant const& t);

// A factory for NavigationFrame objects.
not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const& plugin,
    NavigationFrameParameters const& parameters);

// A factory for PlottingFrame objects.
not_null<std::unique_ptr<PlottingFrame>> NewPlottingFrame(
    Plugin const& plugin,
    PlottingFrameParameters const& parameters);

not_null<OrbitAnalysis*> NewOrbitAnalysis(
    OrbitAnalyser::Analysis* const vessel_analysis,
    Plugin const& plugin,
    int const* const revolutions_per_cycle,
    int const* const days_per_cycle,
    int const ground_track_revolution);

FlightPlan& GetFlightPlan(Plugin const& plugin, char const* vessel_guid);

}  // namespace interface
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
