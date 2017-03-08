#pragma once

#include <limits>

namespace principia {
namespace ksp_plugin {

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
using GUID = std::string;

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = std::uint32_t;

}  // namespace ksp_plugin
}  // namespace principia
