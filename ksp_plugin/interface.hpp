
#pragma once

#include <type_traits>

#include "base/macros.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using base::PullSerializer;
using base::PushDeserializer;
using geometry::Quaternion;
using geometry::R3Element;
using ksp_plugin::NavigationFrame;
using ksp_plugin::Plugin;
using ksp_plugin::RenderedTrajectory;
using ksp_plugin::World;

namespace interface {

struct LineAndIterator {
  explicit LineAndIterator(RenderedTrajectory<World> const& rendered_trajectory)
      : rendered_trajectory(rendered_trajectory) {}
  RenderedTrajectory<World> const rendered_trajectory;
  RenderedTrajectory<World>::const_iterator it;
};

#include "ksp_plugin/interface.generated.h"

extern "C" PRINCIPIA_DLL
void CDECL principia__ActivateRecorder(bool const activate,
                                       bool const verbose);

bool operator==(AdaptiveStepParameters const& left,
                AdaptiveStepParameters const& right);
bool operator==(Burn const& left, Burn const& right);
bool operator==(NavigationFrameParameters const& left,
                NavigationFrameParameters const& right);
bool operator==(NavigationManoeuvre const& left,
                NavigationManoeuvre const& right);
bool operator==(QP const& left, QP const& right);
bool operator==(WXYZ const& left, WXYZ const& right);
bool operator==(XYZ const& left, XYZ const& right);
bool operator==(XYZSegment const& left, XYZSegment const& right);

R3Element<double> ToR3Element(XYZ const& xyz);
WXYZ ToWXYZ(Quaternion const& quaternion);
XYZ ToXYZ(R3Element<double> const& r3_element);

// A factory for NavigationFrame objects.
not_null<std::unique_ptr<NavigationFrame>> NewNavigationFrame(
    Plugin const* const plugin,
    NavigationFrameParameters const& parameters);

}  // namespace interface
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
