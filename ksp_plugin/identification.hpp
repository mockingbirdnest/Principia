#pragma once

#include <limits>

namespace principia {
namespace ksp_plugin {

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = std::uint32_t;

}  // namespace ksp_plugin
}  // namespace principia
