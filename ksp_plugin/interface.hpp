#pragma once

#include <type_traits>

#include "base/macros.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using base::PullSerializer;
using base::PushDeserializer;

namespace ksp_plugin {

struct LineAndIterator {
  explicit LineAndIterator(RenderedTrajectory<World> const& rendered_trajectory)
      : rendered_trajectory(rendered_trajectory) {}
  RenderedTrajectory<World> const rendered_trajectory;
  RenderedTrajectory<World>::const_iterator it;
};

extern "C"
struct XYZ {
  double x, y, z;
};

static_assert(std::is_standard_layout<XYZ>::value,
              "XYZ is used for interfacing");

bool operator==(XYZ const& left, XYZ const& right);

extern "C"
struct XYZSegment {
  XYZ begin, end;
};

static_assert(std::is_standard_layout<XYZSegment>::value,
              "XYZSegment is used for interfacing");

bool operator==(XYZSegment const& left, XYZSegment const& right);

extern "C"
struct WXYZ {
  double w, x, y, z;
};

static_assert(std::is_standard_layout<WXYZ>::value,
              "WXYZ is used for interfacing");

bool operator==(WXYZ const& left, WXYZ const& right);

extern "C"
struct QP {
  XYZ q, p;
};

static_assert(std::is_standard_layout<QP>::value,
              "QP is used for interfacing");

bool operator==(QP const& left, QP const& right);

extern "C"
struct KSPPart {
  // TODO(egg): Y U NO USE QP?
  XYZ world_position;
  XYZ world_velocity;
  double mass;
  XYZ gravitational_acceleration_to_be_applied_by_ksp;
  uint32_t id;
};

static_assert(std::is_standard_layout<KSPPart>::value,
              "KSPPart is used for interfacing");

extern "C"
struct NavigationFrameParameters {
  int extension;  // 6000 to 6999, see serialization::DynamicFrame.
  int centre_index;
  int primary_index;
  int secondary_index;
};

static_assert(std::is_standard_layout<NavigationFrameParameters>::value,
              "NavigationFrameParameters is used for interfacing");

// Sets stderr to log INFO, and redirects stderr, which Unity does not log, to
// "<KSP directory>/stderr.log".  This provides an easily accessible file
// containing a sufficiently verbose log of the latest session, instead of
// requiring users to dig in the archive of all past logs at all severities.
// This archive is written to
// "<KSP directory>/glog/Principia/<SEVERITY>.<date>-<time>.<pid>",
// where date and time are in ISO 8601 basic format.
extern "C" PRINCIPIA_DLL
void CDECL principia__InitGoogleLogging();

// If |activate| is true and there is no active journal, create one and
// activate it.  If |activate| is false and there is an active journal,
// deactivate it.  Does nothing if there is already a journal in the desired
// state.
extern "C" PRINCIPIA_DLL
void CDECL principia__ActivateRecorder(bool const activate);

// Log messages at a level |<= max_severity| are buffered.
// Log messages at a higher level are flushed immediately.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetBufferedLogging(int const max_severity);
extern "C" PRINCIPIA_DLL
int CDECL principia__GetBufferedLogging();
// Sets the maximum number of seconds which logs may be buffered for.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetBufferDuration(int const seconds);
extern "C" PRINCIPIA_DLL
int CDECL principia__GetBufferDuration();
// Log suppression level: messages logged at a lower level than this are
// suppressed.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetSuppressedLogging(int const min_severity);
extern "C" PRINCIPIA_DLL
int CDECL principia__GetSuppressedLogging();
// Show all VLOG(m) messages for |m <= level|.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetVerboseLogging(int const level);
extern "C" PRINCIPIA_DLL
int CDECL principia__GetVerboseLogging();
// Make it so that all log messages of at least |min_severity| are logged to
// stderr (in addition to logging to the usual log file(s)).
extern "C" PRINCIPIA_DLL
void CDECL principia__SetStderrLogging(int const min_severity);
extern "C" PRINCIPIA_DLL
int CDECL principia__GetStderrLogging();

// Exports |LOG(SEVERITY) << message| for fast logging from the C# adapter.
// This will always evaluate its argument even if the corresponding log severity
// is disabled, so it is less efficient than LOG(INFO).  It will not report the
// line and file of the caller.
extern "C" PRINCIPIA_DLL
void CDECL principia__LogInfo(char const* const message);
extern "C" PRINCIPIA_DLL
void CDECL principia__LogWarning(char const* const message);
extern "C" PRINCIPIA_DLL
void CDECL principia__LogError(char const* const message);
extern "C" PRINCIPIA_DLL
void CDECL principia__LogFatal(char const* const message);

// Returns a pointer to a plugin constructed with the arguments given.
// The caller takes ownership of the result.
extern "C" PRINCIPIA_DLL
Plugin* CDECL principia__NewPlugin(
    double const initial_time,
    double const planetarium_rotation_in_degrees);

// Deletes and nulls |*plugin|.
// |plugin| must not be null.  No transfer of ownership of |*plugin|, takes
// ownership of |**plugin|.
extern "C" PRINCIPIA_DLL
void CDECL principia__DeletePlugin(Plugin const** const plugin);

extern "C" PRINCIPIA_DLL
void CDECL principia__DirectlyInsertCelestial(
  Plugin* const plugin,
  int const celestial_index,
  int const* const parent_index,
  char const* const gravitational_parameter,
  char const* const axis_right_ascension,
  char const* const axis_declination,
  char const* const j2,
  char const* const reference_radius,
  char const* const x,
  char const* const y,
  char const* const z,
  char const* const vx,
  char const* const vy,
  char const* const vz);

// Calls |plugin->InsertCelestial| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
void CDECL principia__InsertCelestial(Plugin* const plugin,
                                      int const celestial_index,
                                      double const gravitational_parameter,
                                      int const parent_index,
                                      QP const from_parent);

extern "C" PRINCIPIA_DLL
void CDECL principia__InsertSun(Plugin* const plugin,
                                int const celestial_index,
                                double const gravitational_parameter);

// Calls |plugin->UpdateCelestialHierarchy| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
void CDECL principia__UpdateCelestialHierarchy(Plugin const* const plugin,
                                               int const celestial_index,
                                               int const parent_index);

// Calls |plugin->EndInitialization|.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
void CDECL principia__EndInitialization(Plugin* const plugin);

// Calls |plugin->InsertOrKeepVessel| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
bool CDECL principia__InsertOrKeepVessel(Plugin* const plugin,
                                         char const* const vessel_guid,
                                         int const parent_index);

// Calls |plugin->SetVesselStateOffset| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetVesselStateOffset(Plugin* const plugin,
                                           char const* const vessel_guid,
                                           QP const from_parent);

extern "C" PRINCIPIA_DLL
void CDECL principia__AdvanceTime(Plugin* const plugin,
                                  double const t,
                                  double const planetarium_rotation);

extern "C" PRINCIPIA_DLL
void CDECL principia__ForgetAllHistoriesBefore(Plugin* const plugin,
                                               double const t);

// Calls |plugin->VesselFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
QP CDECL principia__VesselFromParent(Plugin const* const plugin,
                                     char const* const vessel_guid);

// Calls |plugin->CelestialFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
QP CDECL principia__CelestialFromParent(Plugin const* const plugin,
                                        int const celestial_index);

// Calls |plugin->NewBodyCentredNonRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
// TODO(phl): The parameter should be named |centre_index|.
extern "C" PRINCIPIA_DLL
NavigationFrame* CDECL principia__NewBodyCentredNonRotatingNavigationFrame(
    Plugin const* const plugin,
    int const reference_body_index);

// Calls |plugin->NewBarycentricRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
extern "C" PRINCIPIA_DLL
NavigationFrame* CDECL principia__NewBarycentricRotatingNavigationFrame(
    Plugin const* const plugin,
    int const primary_index,
    int const secondary_index);

// Calls |plugin| to create a |NavigationFrame| using the given |parameters|.
extern "C" PRINCIPIA_DLL
NavigationFrame* CDECL principia__NewNavigationFrame(
    Plugin const* const plugin,
    NavigationFrameParameters const parameters);

// Returns the parameters for the given frame, which must not be null.
extern "C" PRINCIPIA_DLL
NavigationFrameParameters CDECL principia__GetNavigationFrameParameters(
    NavigationFrame const* const navigation_frame);

// |navigation_frame| must not be null.  No transfer of ownership of
// |*navigation_frame|, takes ownership of |**navigation_frame|, nulls
// |*navigation_frame|.
extern "C" PRINCIPIA_DLL
void CDECL principia__SetPlottingFrame(
    Plugin* const plugin,
    NavigationFrame** const navigation_frame);

// Returns the frame last set by |plugin->SetPlottingFrame|.  No transfer of
// ownership.  The returned pointer is never null.
extern "C" PRINCIPIA_DLL
NavigationFrame const* CDECL principia__GetPlottingFrame(
    Plugin const* const plugin);

extern "C" PRINCIPIA_DLL
void principia__UpdatePrediction(Plugin const* const plugin,
                                 char const* const vessel_guid);

// Returns the result of |plugin->RenderedVesselTrajectory| called with the
// arguments given, together with an iterator to its beginning.
// |plugin| must not be null.  No transfer of ownership of |plugin|.  The caller
// gets ownership of the result.  |frame| must not be null.  No transfer of
// ownership of |frame|.
extern "C" PRINCIPIA_DLL
LineAndIterator* CDECL principia__RenderedVesselTrajectory(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position);

extern "C" PRINCIPIA_DLL
bool principia__HasPrediction(Plugin const* const plugin,
                              char const* const vessel_guid);

extern "C" PRINCIPIA_DLL
LineAndIterator* CDECL principia__RenderedPrediction(
    Plugin* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position);

extern "C" PRINCIPIA_DLL
int principia__FlightPlanSize(Plugin const* const plugin,
                              char const* const vessel_guid);

extern "C" PRINCIPIA_DLL
LineAndIterator* CDECL principia__RenderedFlightPlan(
    Plugin* const plugin,
    char const* const vessel_guid,
    int const plan_phase,
    XYZ const sun_world_position);

extern "C" PRINCIPIA_DLL
void CDECL principia__SetPredictionLength(Plugin* const plugin,
                                          double const t);

extern "C" PRINCIPIA_DLL
void CDECL principia__SetPredictionLengthTolerance(Plugin* const plugin,
                                                   double const l);

extern "C" PRINCIPIA_DLL
void CDECL principia__SetPredictionSpeedTolerance(Plugin* const plugin,
                                                  double const v);

extern "C" PRINCIPIA_DLL
bool CDECL principia__HasVessel(Plugin* const plugin,
                                char const* const vessel_guid);

// Returns |line_and_iterator->rendered_trajectory.size()|.
// |line_and_iterator| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
int CDECL principia__NumberOfSegments(
    LineAndIterator const* const line_and_iterator);

// Returns the |XYZSegment| corresponding to the |LineSegment|
// |*line_and_iterator->it|, then increments |line_and_iterator->it|.
// |line_and_iterator| must not be null.  |line_and_iterator->it| must not be
// the end of |line_and_iterator->rendered_trajectory|.  No transfer of
// ownership.
extern "C" PRINCIPIA_DLL
XYZSegment CDECL principia__FetchAndIncrement(
    LineAndIterator* const line_and_iterator);

// Returns |true| if and only if |line_and_iterator->it| is the end of
// |line_and_iterator->rendered_trajectory|.
// |line_and_iterator| must not be null.  No transfer of ownership.
extern "C" PRINCIPIA_DLL
bool CDECL principia__AtEnd(LineAndIterator const* const line_and_iterator);

// Deletes and nulls |*line_and_iterator|.
// |line_and_iterator| must not be null.  No transfer of ownership of
// |*line_and_iterator|, takes ownership of |**line_and_iterator|.
extern "C" PRINCIPIA_DLL
void CDECL principia__DeleteLineAndIterator(
    LineAndIterator** const line_and_iterator);

extern "C" PRINCIPIA_DLL
void CDECL principia__AddVesselToNextPhysicsBubble(
    Plugin* const plugin,
    char const* const vessel_guid,
    KSPPart const* const parts,
    int count);

extern "C" PRINCIPIA_DLL
bool CDECL principia__PhysicsBubbleIsEmpty(Plugin const* const plugin);

extern "C" PRINCIPIA_DLL
XYZ CDECL principia__BubbleDisplacementCorrection(Plugin const* const plugin,
                                                  XYZ const sun_position);

extern "C" PRINCIPIA_DLL
XYZ CDECL principia__BubbleVelocityCorrection(Plugin const* const plugin,
                                              int const reference_body_index);

extern "C" PRINCIPIA_DLL
WXYZ CDECL principia__NavballOrientation(
    Plugin const* const plugin,
    XYZ const sun_world_position,
    XYZ const ship_world_position);

extern "C" PRINCIPIA_DLL
XYZ CDECL principia__VesselTangent(Plugin const* const plugin,
                                   char const* const vessel_guid);

extern "C" PRINCIPIA_DLL
XYZ CDECL principia__VesselNormal(Plugin const* const plugin,
                                  char const* const vessel_guid);

extern "C" PRINCIPIA_DLL
XYZ CDECL principia__VesselBinormal(Plugin const* const plugin,
                                    char const* const vessel_guid);

extern "C" PRINCIPIA_DLL
double CDECL principia__CurrentTime(Plugin const* const plugin);

// |plugin| must not be null.  The caller takes ownership of the result, except
// when it is null (at the end of the stream).  No transfer of ownership of
// |*plugin|.  |*serializer| must be null on the first call and must be passed
// unchanged to the successive calls; its ownership is not transferred.
extern "C" PRINCIPIA_DLL
char const* CDECL principia__SerializePlugin(
    Plugin const* const plugin,
    PullSerializer** const serializer);

// Deletes and nulls |*serialization|.
// |serialization| must not be null.  No transfer of ownership of
// |*serialization|, takes ownership of |**serialization|.
extern "C" PRINCIPIA_DLL
void CDECL principia__DeletePluginSerialization(
    char const** const serialization);

// The caller takes ownership of |**plugin| when it is not null.  No transfer of
// ownership of |*serialization| or |**deserializer|.  |*deserializer| and
// |*plugin| must be null on the first call and must be passed unchanged to the
// successive calls.  The caller must perform an extra call with
// |serialization_size| set to 0 to indicate the end of the input stream.  When
// this last call returns, |*plugin| is not null and may be used by the caller.
extern "C" PRINCIPIA_DLL
void CDECL principia__DeserializePlugin(
    char const* const serialization,
    int const serialization_size,
    PushDeserializer** const deserializer,
    Plugin const** const plugin);

// Says hello, convenient for checking that calls to the DLL work.
extern "C" PRINCIPIA_DLL
char const* CDECL principia__SayHello();

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
