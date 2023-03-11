#pragma once

#include <string>
#include <typeindex>
#include <type_traits>
#include <utility>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "base/macros.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/iterators.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/planetarium.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace interface {

// This is used for interfacing, and should only appear in C++ code in tests
// and generated code; we allow ourselves to pollute the |interface| namespace
// with convenience |using|s.

using namespace principia::base::_not_null;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_iterators;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_dynamic_frame;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Takes ownership of |**pointer| and returns it to the caller.  Nulls
// |*pointer|.  |pointer| must not be null.  No transfer of ownership of
// |*pointer|.
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
bool operator==(OrbitAnalysis const& left, OrbitAnalysis const& right);
bool operator==(OrbitGroundTrack const& left, OrbitGroundTrack const& right);
bool operator==(OrbitRecurrence const& left, OrbitRecurrence const& right);
bool operator==(OrbitalElements const& left, OrbitalElements const& right);
bool operator==(QP const& left, QP const& right);
bool operator==(QPRW const& left, QPRW const& right);
bool operator==(WXYZ const& left, WXYZ const& right);
bool operator==(XY const& left, XY const& right);
bool operator==(XYZ const& left, XYZ const& right);

// Conversions between interchange data and typed data.

physics::Ephemeris<Barycentric>::AdaptiveStepParameters
FromAdaptiveStepParameters(
    AdaptiveStepParameters const& adaptive_step_parameters);
std::pair<physics::Ephemeris<Barycentric>::AdaptiveStepParameters,
          physics::Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters>
FromFlightPlanAdaptiveStepParameters(FlightPlanAdaptiveStepParameters const&
                                         flight_plan_adaptive_step_parameters);

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
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters);
FlightPlanAdaptiveStepParameters ToFlightPlanAdaptiveStepParameters(
    physics::Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    physics::Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters);

KeplerianElements ToKeplerianElements(
    physics::KeplerianElements<Barycentric> const& keplerian_elements);

QP ToQP(DegreesOfFreedom<World> const& dof);
QP ToQP(RelativeDegreesOfFreedom<AliceSun> const& relative_dof);

// Ownership of the status and its message is transferred to the caller.
Status* ToNewStatus(absl::Status const& status);

WXYZ ToWXYZ(geometry::Quaternion const& quaternion);

XY ToXY(geometry::RP2Point<Length, Camera> const& rp2_point);

XYZ ToXYZ(geometry::R3Element<double> const& r3_element);
XYZ ToXYZ(Position<World> const& position);
XYZ ToXYZ(Vector<double, World> const& direction);
XYZ ToXYZ(Velocity<Frenet<NavigationFrame>> const& velocity);
XYZ ToXYZ(Bivector<AngularMomentum, World> const& angular_momentum);

// Conversions between interchange data and typed data that depend on the state
// of the plugin.
Instant FromGameTime(Plugin const& plugin, double t);
double ToGameTime(Plugin const& plugin, Instant const& t);

// A factory for NavigationFrame objects.
not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const& plugin,
    NavigationFrameParameters const& parameters);

}  // namespace interface
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
