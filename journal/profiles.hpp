#pragma once

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "base/push_pull_callback.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"  // 🧙 For symbols in interface.
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/planetarium.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/quantities.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

// This file is not included from other headers, only translation units, so we
// allow pollution of the `journal` namespace.
using interface::AdaptiveStepParameters;
using interface::BodyGeopotentialElement;
using interface::BodyParameters;
using interface::Burn;
using interface::ConfigurationAccuracyParameters;
using interface::ConfigurationAdaptiveStepParameters;
using interface::ConfigurationDownsamplingParameters;
using interface::ConfigurationFixedStepParameters;
using interface::EquatorialCrossings;
using interface::FlightPlanAdaptiveStepParameters;
using interface::Interval;
using interface::Iterator;
using interface::KeplerianElements;
using interface::NavigationFrameParameters;
using interface::NavigationManoeuvre;
using interface::NavigationManoeuvreFrenetTrihedron;
using interface::Node;
using interface::OrbitalElements;
using interface::OrbitAnalysis;
using interface::OrbitRecurrence;
using interface::Origin;
using interface::PlottingFrameParameters;
using interface::QP;
using interface::QPRW;
using interface::SolarTimesOfNodes;
using interface::Status;
using interface::TQP;
using interface::WXYZ;
using interface::XY;
using interface::XYZ;
using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;
using namespace principia::base::_push_pull_callback;
using namespace principia::journal::_player;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_vessel;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_quantities;

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
