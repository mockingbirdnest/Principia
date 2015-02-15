#pragma once

#include <type_traits>

#include "base/macros.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/transforms.hpp"

namespace principia {
namespace ksp_plugin {

struct LineAndIterator {
  explicit LineAndIterator(RenderedTrajectory<World> const& rendered_trajectory)
      : rendered_trajectory(rendered_trajectory) {}
  RenderedTrajectory<World> const rendered_trajectory;
  RenderedTrajectory<World>::const_iterator it;
};

}  // namespace ksp_plugin
}  // namespace principia

using principia::ksp_plugin::Barycentric;
using principia::ksp_plugin::LineAndIterator;
using principia::ksp_plugin::Plugin;
using principia::ksp_plugin::Rendering;
using principia::physics::Transforms;

extern "C"
struct XYZ {
  double x, y, z;
};

static_assert(std::is_standard_layout<XYZ>::value,
              "XYZ is used for interfacing");

extern "C"
struct XYZSegment {
  XYZ begin, end;
};

static_assert(std::is_standard_layout<XYZSegment>::value,
              "XYZSegment is used for interfacing");

extern "C"
struct QP {
  XYZ q, p;
};

static_assert(std::is_standard_layout<QP>::value,
              "QP is used for interfacing");

extern "C"
struct KSPPart {
  XYZ world_position;
  XYZ world_velocity;
  double mass;
  XYZ gravitational_acceleration_to_be_applied_by_ksp;
  uint32_t id;
};

static_assert(std::is_standard_layout<KSPPart>::value,
              "KSPPart is used for interfacing");

// Sets stderr to log INFO, and redirects stderr, which Unity does not log, to
// "<KSP directory>/stderr.log".  This provides an easily accessible file
// containing a sufficiently verbose log of the latest session, instead of
// requiring users to dig in the archive of all past logs at all severities.
// This archive is written to
// "<KSP directory>/glog/Principia/<SEVERITY>.<date>-<time>.<pid>",
// where date and time are in ISO 8601 basic format.
extern "C" DLLEXPORT
void CDECL principia__InitGoogleLogging();

// Log messages at a level |<= max_severity| are buffered.
// Log messages at a higher level are flushed immediately.
extern "C" DLLEXPORT
void CDECL principia__SetBufferedLogging(int const max_severity);
extern "C" DLLEXPORT
int CDECL principia__GetBufferedLogging();
// Sets the maximum number of seconds which logs may be buffered for.
extern "C" DLLEXPORT
void CDECL principia__SetBufferDuration(int const seconds);
extern "C" DLLEXPORT
int CDECL principia__GetBufferDuration();
// Log suppression level: messages logged at a lower level than this are
// suppressed.
extern "C" DLLEXPORT
void CDECL principia__SetSuppressedLogging(int const min_severity);
extern "C" DLLEXPORT
int CDECL principia__GetSuppressedLogging();
// Show all VLOG(m) messages for |m <= level|.
extern "C" DLLEXPORT
void CDECL principia__SetVerboseLogging(int const level);
extern "C" DLLEXPORT
int CDECL principia__GetVerboseLogging();
// Make it so that all log messages of at least |min_severity| are logged to
// stderr (in addition to logging to the usual log file(s)).
extern "C" DLLEXPORT
void CDECL principia__SetStderrLogging(int const min_severity);
extern "C" DLLEXPORT
int CDECL principia__GetStderrLogging();

// Exports |LOG(SEVERITY) << message| for fast logging from the C# adapter.
// This will always evaluate its argument even if the corresponding log severity
// is disabled, so it is less efficient than LOG(INFO).  It will not report the
// line and file of the caller.
extern "C" DLLEXPORT
void CDECL principia__LogInfo(char const* message);
extern "C" DLLEXPORT
void CDECL principia__LogWarning(char const* message);
extern "C" DLLEXPORT
void CDECL principia__LogError(char const* message);
extern "C" DLLEXPORT
void CDECL principia__LogFatal(char const* message);

// Returns a pointer to a plugin constructed with the arguments given.
// The caller takes ownership of the result.
extern "C" DLLEXPORT
Plugin* CDECL principia__NewPlugin(
    double const initial_time,
    int const sun_index,
    double const sun_gravitational_parameter,
    double const planetarium_rotation_in_degrees);

// Deletes and nulls |*plugin|.
// |plugin| must not be null.  No transfer of ownership of |*plugin|, takes
// ownership of |**plugin|.
extern "C" DLLEXPORT
void CDECL principia__DeletePlugin(Plugin const** const plugin);

// Calls |plugin->InsertCelestial| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
void CDECL principia__InsertCelestial(Plugin* const plugin,
                                      int const celestial_index,
                                      double const gravitational_parameter,
                                      int const parent_index,
                                      QP const from_parent);

// Calls |plugin->UpdateCelestialHierarchy| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
void CDECL principia__UpdateCelestialHierarchy(Plugin const* const plugin,
                                               int const celestial_index,
                                               int const parent_index);

// Calls |plugin->EndInitialization|.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
void CDECL principia__EndInitialization(Plugin* const plugin);

// Calls |plugin->InsertOrKeepVessel| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
bool CDECL principia__InsertOrKeepVessel(Plugin* const plugin,
                                         char const* vessel_guid,
                                         int const parent_index);

// Calls |plugin->SetVesselStateOffset| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
void CDECL principia__SetVesselStateOffset(Plugin* const plugin,
                                           char const* vessel_guid,
                                           QP const from_parent);

extern "C" DLLEXPORT
void CDECL principia__AdvanceTime(Plugin* const plugin,
                                  double const t,
                                  double const planetarium_rotation);

// Calls |plugin->VesselFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
QP CDECL principia__VesselFromParent(Plugin const* const plugin,
                                     char const* vessel_guid);

// Calls |plugin->CelestialFromParent| with the arguments given.
// |plugin| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
QP CDECL principia__CelestialFromParent(Plugin const* const plugin,
                                        int const celestial_index);

// Calls |plugin->NewBodyCentredNonRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
extern "C" DLLEXPORT
Transforms<Barycentric, Rendering, Barycentric>* CDECL
principia__NewBodyCentredNonRotatingTransforms(Plugin const* const plugin,
                                               int const reference_body_index);

// Calls |plugin->NewBarycentricRotatingFrame| with the arguments given.
// |plugin| must not be null.  The caller gets ownership of the returned object.
extern "C" DLLEXPORT
Transforms<Barycentric, Rendering, Barycentric>* CDECL
principia__NewBarycentricRotatingTransforms(Plugin const* const plugin,
                                            int const primary_index,
                                            int const secondary_index);

// Deletes and nulls |*transforms|.
// |transforms| must not be null.  No transfer of ownership of |*transforms|,
// takes ownership of |**transforms|.
extern "C" DLLEXPORT
void CDECL principia__DeleteTransforms(
    Transforms<Barycentric, Rendering, Barycentric>** const transforms);

// Returns the result of |plugin->RenderedVesselTrajectory| called with the
// arguments given, together with an iterator to its beginning.
// |plugin| must not be null.  No transfer of ownership of |plugin|.  The caller
// gets ownership of the result.  |frame| must not be null.  No transfer of
// ownership of |frame|.
extern "C" DLLEXPORT
LineAndIterator* CDECL principia__RenderedVesselTrajectory(
    Plugin const* const plugin,
    char const* vessel_guid,
    Transforms<Barycentric, Rendering, Barycentric>* const transforms,
    XYZ const sun_world_position);

// Returns |line_and_iterator->rendered_trajectory.size()|.
// |line_and_iterator| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
int CDECL principia__NumberOfSegments(LineAndIterator const* line_and_iterator);

// Returns the |XYZSegment| corresponding to the |LineSegment|
// |*line_and_iterator->it|, then increments |line_and_iterator->it|.
// |line_and_iterator| must not be null.  |line_and_iterator->it| must not be
// the end of |line_and_iterator->rendered_trajectory|.  No transfer of
// ownership.
extern "C" DLLEXPORT
XYZSegment CDECL principia__FetchAndIncrement(
    LineAndIterator* const line_and_iterator);

// Returns |true| if and only if |line_and_iterator->it| is the end of
// |line_and_iterator->rendered_trajectory|.
// |line_and_iterator| must not be null.  No transfer of ownership.
extern "C" DLLEXPORT
bool CDECL principia__AtEnd(LineAndIterator* const line_and_iterator);

// Deletes and nulls |*line_and_iterator|.
// |line_and_iterator| must not be null.  No transfer of ownership of
// |*line_and_iterator|, takes ownership of |**line_and_iterator|.
extern "C" DLLEXPORT
void CDECL principia__DeleteLineAndIterator(
    LineAndIterator** const line_and_iterator);

extern "C" DLLEXPORT
XYZ CDECL principia__VesselWorldPosition(Plugin const* const plugin,
                                         char const* vessel_guid,
                                         XYZ const parent_world_position);

extern "C" DLLEXPORT
XYZ CDECL principia__VesselWorldVelocity(Plugin const* const plugin,
                                         char const* vessel_guid,
                                         XYZ const parent_world_velocity,
                                         double const parent_rotation_period);

extern "C" DLLEXPORT
void CDECL principia__AddVesselToNextPhysicsBubble(Plugin* const plugin,
                                                   char const* vessel_guid,
                                                   KSPPart const* const parts,
                                                   int count);

extern "C" DLLEXPORT
bool CDECL principia__PhysicsBubbleIsEmpty(Plugin const* const plugin);

extern "C" DLLEXPORT
XYZ CDECL principia__BubbleDisplacementCorrection(Plugin const* const plugin,
                                                  XYZ const sun_position);

extern "C" DLLEXPORT
XYZ CDECL principia__BubbleVelocityCorrection(Plugin const* const plugin,
                                              int const reference_body_index);

extern "C" DLLEXPORT
double CDECL principia__current_time(Plugin const* const plugin);

// |plugin| must not be null.  The caller takes ownership of the result.  No
// transfer of ownership of |*plugin|.
extern "C" DLLEXPORT
char const* CDECL principia__SerializePlugin(Plugin const* const plugin);

// Deletes and nulls |*serialization|.
// |serialization| must not be null.  No transfer of ownership of
// |*serialization|, takes ownership of |**serialization|.
extern "C" DLLEXPORT
void CDECL principia__DeletePluginSerialization(
    char const** const serialization);

// The caller takes ownership of the result.  No transfer of ownership of
// |*serialization|.
extern "C" DLLEXPORT
Plugin* CDECL principia__DeserializePlugin(char const* const serialization,
                                           int const serialization_size);

// Says hello, convenient for checking that calls to the DLL work.
extern "C" DLLEXPORT
char const* CDECL principia__SayHello();

#undef CDECL
#undef DLLEXPORT
