#pragma once

#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;
using interface::KSPPart;
using interface::LineAndIterator;
using interface::NavigationFrameParameters;
using interface::QP;
using interface::WXYZ;
using interface::XYZ;
using interface::XYZSegment;

namespace journal {

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
