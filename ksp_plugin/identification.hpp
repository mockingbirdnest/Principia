#pragma once

#include <limits>

namespace principia {
namespace ksp_plugin {

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = std::uint32_t;

// An ID given to dummy parts used by vessels that appear unloaded.  Note that
// nothing prevents an actual part from having this ID.
constexpr PartId DummyPartId = std::numeric_limits<PartId>::max();

}  // namespace ksp_plugin
}  // namespace principia
