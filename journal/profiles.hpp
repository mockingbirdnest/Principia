#pragma once

#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::not_null;
using ksp_plugin::KSPPart;
using ksp_plugin::LineAndIterator;
using ksp_plugin::NavigationFrame;
using ksp_plugin::Plugin;
using ksp_plugin::QP;
using ksp_plugin::WXYZ;
using ksp_plugin::XYZ;
using ksp_plugin::XYZSegment;

namespace journal {

#include "journal/profiles.gen.hpp"

}  // namespace journal
}  // namespace principia
