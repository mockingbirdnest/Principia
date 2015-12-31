#include "journal/profiles.hpp"

#include <fstream>
#include <list>
#include <string>
#include <type_traits>
#include <vector>

#include "base/map_util.hpp"
#include "glog/logging.h"

namespace principia {

using base::FindOrDie;

namespace journal {
namespace {

template<typename T>
void Insert(not_null<Player::PointerMap*> const pointer_map,
            std::uint64_t const address,
            T* const pointer) {
  void* const inserted_pointer =
      static_cast<void*>(const_cast<std::remove_cv<T>::type*>(pointer));
  auto inserted = pointer_map->emplace(address, inserted_pointer);
  if (!inserted.second) {
    CHECK_EQ(inserted.first->second, inserted_pointer);
  }
}

void Delete(not_null<Player::PointerMap*> const pointer_map,
            std::uint64_t const address) {
  if (reinterpret_cast<void*>(address) != nullptr) {
    auto const it = pointer_map->find(address);
    CHECK(it != pointer_map->end()) << address;
    pointer_map->erase(it);
  }
}

template<typename T,
         typename = typename std::enable_if<std::is_pointer<T>::value>::type>
T DeserializePointer(Player::PointerMap const& pointer_map,
                     std::uint64_t const address) {
  if (reinterpret_cast<T>(address) == nullptr) {
    return nullptr;
  } else {
    return reinterpret_cast<T>(FindOrDie(pointer_map, address));
  }
}

WXYZ DeserializeWXYZ(serialization::WXYZ const& wxyz) {
  return {wxyz.w(), wxyz.x(), wxyz.y(), wxyz.z()};
}

XYZ DeserializeXYZ(serialization::XYZ const& xyz) {
  return {xyz.x(), xyz.y(), xyz.z()};
}

XYZSegment DeserializeXYZSegment(serialization::XYZSegment const& xyz_segment) {
  return {DeserializeXYZ(xyz_segment.begin()),
          DeserializeXYZ(xyz_segment.end())};
}

QP DeserializeQP(serialization::QP const& qp) {
  return {DeserializeXYZ(qp.q()), DeserializeXYZ(qp.p())};
}

KSPPart DeserializeKSPPart(serialization::KSPPart const& ksp_part) {
  return {DeserializeXYZ(ksp_part.world_position()),
          DeserializeXYZ(ksp_part.world_velocity()),
          ksp_part.mass(),
          DeserializeXYZ(
              ksp_part.gravitational_acceleration_to_be_applied_by_ksp()),
          ksp_part.id()};
}

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

serialization::WXYZ SerializeWXYZ(WXYZ const& wxyz) {
  serialization::WXYZ m;
  m.set_w(wxyz.w);
  m.set_x(wxyz.x);
  m.set_y(wxyz.y);
  m.set_z(wxyz.z);
  return m;
}

serialization::XYZ SerializeXYZ(XYZ const& xyz) {
  serialization::XYZ m;
  m.set_x(xyz.x);
  m.set_y(xyz.y);
  m.set_z(xyz.z);
  return m;
}

serialization::XYZSegment SerializeXYZSegment(XYZSegment const& xyz_segment) {
  serialization::XYZSegment m;
  *m.mutable_begin() = SerializeXYZ(xyz_segment.begin);
  *m.mutable_end() = SerializeXYZ(xyz_segment.end);
  return m;
}

serialization::QP SerializeQP(QP const& qp) {
  serialization::QP m;
  *m.mutable_p() = SerializeXYZ(qp.p);
  *m.mutable_q() = SerializeXYZ(qp.q);
  return m;
}

serialization::KSPPart SerializeKSPPart(KSPPart const& ksp_part) {
  serialization::KSPPart m;
  *m.mutable_world_position() = SerializeXYZ(ksp_part.world_position);
  *m.mutable_world_velocity() = SerializeXYZ(ksp_part.world_velocity);
  m.set_mass(ksp_part.mass);
  *m.mutable_gravitational_acceleration_to_be_applied_by_ksp() =
      SerializeXYZ(ksp_part.gravitational_acceleration_to_be_applied_by_ksp);
  m.set_id(ksp_part.id);
  return m;
}

}  // namespace

#include "journal/profiles.gen.cpp"

}  // namespace journal
}  // namespace principia
