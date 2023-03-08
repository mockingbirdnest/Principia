#include "ksp_plugin/identification.hpp"

#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace _identification {
namespace internal {

bool PartByPartIdComparator::operator()(not_null<Part*> const left,
                                        not_null<Part*> const right) const {
  return left->part_id() < right->part_id();
}

bool PartByPartIdComparator::operator()(
    not_null<Part const*> const left,
    not_null<Part const*> const right) const {
  return left->part_id() < right->part_id();
}

bool VesselByGUIDComparator::operator()(not_null<Vessel*> const left,
                                        not_null<Vessel*> const right) const {
  return left->guid() < right->guid();
}

bool VesselByGUIDComparator::operator()(
    not_null<Vessel const*> const left,
    not_null<Vessel const*> const right) const {
  return left->guid() < right->guid();
}

}  // namespace internal
}  // namespace _identification
}  // namespace ksp_plugin
}  // namespace principia
