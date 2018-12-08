
#pragma once

#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

// This file is not included from other headers, only translation units, so we
// allow pollution of the principia namespace.
using base::not_constructible;
using base::not_null;
using base::PullSerializer;
using base::PushDeserializer;
using interface::AdaptiveStepParameters;
using interface::BodyGeopotentialElement;
using interface::BodyParameters;
using interface::Burn;
using interface::ConfigurationAccuracyParameters;
using interface::ConfigurationAdaptiveStepParameters;
using interface::ConfigurationFixedStepParameters;
using interface::FlightPlanAdaptiveStepParameters;
using interface::KeplerianElements;
using interface::Iterator;
using interface::NavigationFrameParameters;
using interface::NavigationManoeuvre;
using interface::NavigationManoeuvreFrenetTrihedron;
using interface::Origin;
using interface::QP;
using interface::Status;
using interface::WXYZ;
using interface::XY;
using interface::XYZ;
using ksp_plugin::NavigationFrame;
using ksp_plugin::PileUpFuture;
using ksp_plugin::Planetarium;
using ksp_plugin::Plugin;
using ksp_plugin::Vessel;

namespace journal {

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
