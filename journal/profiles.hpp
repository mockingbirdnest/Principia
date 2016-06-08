
#pragma once

#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;
using interface::AdaptiveStepParameters;
using interface::Burn;
using interface::KeplerianElements;
using interface::Iterator;
using interface::KSPPart;
using interface::NavigationFrameParameters;
using interface::NavigationManoeuvre;
using interface::QP;
using interface::WXYZ;
using interface::XYZ;

namespace journal {

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
