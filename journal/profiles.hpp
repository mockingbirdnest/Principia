#pragma once

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

// This file is not included from other headers, only translation units, so we
// allow pollution of the |journal| namespace.
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
using interface::OrbitalElements;
using interface::OrbitAnalysis;
using interface::OrbitRecurrence;
using interface::Origin;
using interface::QP;
using interface::QPRW;
using interface::SolarTimesOfNodes;
using interface::Status;
using interface::WXYZ;
using interface::XY;
using interface::XYZ;
using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_pull_serializer;
using namespace principia::base::_push_deserializer;
using namespace principia::journal::_player;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::ksp_plugin::_plugin;
using namespace principia::ksp_plugin::_vessel;

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
