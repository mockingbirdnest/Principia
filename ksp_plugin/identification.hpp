#pragma once

#include <limits>
#include <map>
#include <set>
#include <string>

#include "base/not_null.hpp"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE_FROM(part, class, Part);
FORWARD_DECLARE_FROM(vessel, class, Vessel);

namespace internal_identification {

using namespace principia::base::_not_null;

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
using GUID = std::string;

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = std::uint32_t;

// Comparator by PartId.  Useful for ensuring a consistent ordering in sets of
// pointers to Parts.
struct PartByPartIdComparator {
  bool operator()(not_null<Part*> left, not_null<Part*> right) const;
  bool operator()(not_null<Part const*> left,
                  not_null<Part const*> right) const;
};

// Comparator by GUID.  Useful for ensuring a consistent ordering in sets of
// pointers to Vessels.
struct VesselByGUIDComparator {
  bool operator()(not_null<Vessel*> left, not_null<Vessel*> right) const;
  bool operator()(not_null<Vessel const*> left,
                  not_null<Vessel const*> right) const;
};

template<typename T>
using PartTo = std::map<not_null<Part*>,
                        T,
                        PartByPartIdComparator>;
using VesselSet = std::set<not_null<Vessel*>,
                           VesselByGUIDComparator>;
using VesselConstSet = std::set<not_null<Vessel const*>,
                                VesselByGUIDComparator>;

}  // namespace internal_identification

using internal_identification::GUID;
using internal_identification::PartByPartIdComparator;
using internal_identification::PartId;
using internal_identification::PartTo;
using internal_identification::VesselByGUIDComparator;
using internal_identification::VesselConstSet;
using internal_identification::VesselSet;

}  // namespace ksp_plugin
}  // namespace principia
