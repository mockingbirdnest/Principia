#pragma once

#include <type_traits>

#include "base/macros.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {

using base::PullSerializer;
using base::PushDeserializer;
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
void CDECL principia__ActivateRecorder(bool const activate);

bool operator==(NavigationFrameParameters const& left,
                NavigationFrameParameters const& right);
bool operator==(QP const& left, QP const& right);
bool operator==(WXYZ const& left, WXYZ const& right);
bool operator==(XYZ const& left, XYZ const& right);
bool operator==(XYZSegment const& left, XYZSegment const& right);

}  // namespace interface
}  // namespace principia

#include "ksp_plugin/interface_body.hpp"
