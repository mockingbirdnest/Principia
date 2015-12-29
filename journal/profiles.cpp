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

NavigationFrameParameters DeserializeNavigationFrameParameters(
    serialization::NavigationFrameParameters const& parameters) {
  return {parameters.extension(),
          parameters.centre_index(),
          parameters.primary_index(),
          parameters.secondary_index()};
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

serialization::NavigationFrameParameters SerializeNavigationFrameParameters(
    NavigationFrameParameters const& parameters) {
  serialization::NavigationFrameParameters m;
  m.set_extension(parameters.extension);
  m.set_centre_index(parameters.centre_index);
  m.set_primary_index(parameters.primary_index);
  m.set_secondary_index(parameters.secondary_index);
  return m;
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

void InitGoogleLogging::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {}

void SetBufferedLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_max_severity(in.max_severity);
}

void SetBufferedLogging::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__SetBufferedLogging(message.in().max_severity());
}

void GetBufferedLogging::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_get_buffered_logging(result);
}

void GetBufferedLogging::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_buffered_logging(),
           ksp_plugin::principia__GetBufferedLogging());
}

void SetBufferDuration::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_seconds(in.seconds);
}

void SetBufferDuration::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__SetBufferDuration(message.in().seconds());
}

void GetBufferDuration::Fill(Return const& result,
                             not_null<Message*> const message) {
  message->mutable_return_()->set_get_buffer_duration(result);
}

void GetBufferDuration::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_buffer_duration(),
           ksp_plugin::principia__GetBufferDuration());
}

void SetSuppressedLogging::Fill(In const& in,
                                not_null<Message*> const message) {
  message->mutable_in()->set_min_severity(in.min_severity);
}

void SetSuppressedLogging::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__SetSuppressedLogging(message.in().min_severity());
}

void GetSuppressedLogging::Fill(Return const& result,
                                not_null<Message*> const message) {
  message->mutable_return_()->set_get_suppressed_logging(result);
}

void GetSuppressedLogging::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_suppressed_logging(),
           ksp_plugin::principia__GetSuppressedLogging());
}

void SetVerboseLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_level(in.level);
}

void SetVerboseLogging::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__SetVerboseLogging(message.in().level());
}

void GetVerboseLogging::Fill(Return const& result,
                             not_null<Message*> const message) {
  message->mutable_return_()->set_get_verbose_logging(result);
}

void GetVerboseLogging::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_verbose_logging(),
           ksp_plugin::principia__GetVerboseLogging());
}

void SetStderrLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_min_severity(in.min_severity);
}

void SetStderrLogging::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__SetStderrLogging(message.in().min_severity());
}

void GetStderrLogging::Fill(Return const& result,
                            not_null<Message*> const message) {
  message->mutable_return_()->set_get_stderr_logging(result);
}

void GetStderrLogging::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_stderr_logging(),
           ksp_plugin::principia__GetStderrLogging());
}

void LogInfo::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogInfo::Run(Message const& message,
                  not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__LogInfo(message.in().message().c_str());
}

void LogWarning::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogWarning::Run(Message const& message,
                     not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__LogWarning(message.in().message().c_str());
}

void LogError::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogError::Run(Message const& message,
                   not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__LogError(message.in().message().c_str());
}

void LogFatal::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogFatal::Run(Message const& message,
                   not_null<Player::PointerMap*> const pointer_map) {
  ksp_plugin::principia__LogFatal(message.in().message().c_str());
}

void NewPlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* mutable_in = message->mutable_in();
  mutable_in->set_initial_time(in.initial_time);
  mutable_in->set_planetarium_rotation_in_degrees(
      in.planetarium_rotation_in_degrees);
}

void NewPlugin::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_new_plugin(SerializePointer(result));
}

void NewPlugin::Run(Message const& message,
                    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = ksp_plugin::principia__NewPlugin(
                     in.initial_time(),
                     in.planetarium_rotation_in_degrees());
  Insert(pointer_map, message.return_().new_plugin(), plugin);
}

void DeletePlugin::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_plugin(SerializePointer(*in.plugin));
}

void DeletePlugin::Fill(Out const& out, not_null<Message*> const message) {
  message->mutable_out()->set_plugin(SerializePointer(*out.plugin));
}

void DeletePlugin::Run(Message const& message,
                       not_null<Player::PointerMap*> const pointer_map) {
  auto* plugin = DeserializePointer<Plugin const*>(*pointer_map,
                                                   message.in().plugin());
  ksp_plugin::principia__DeletePlugin(&plugin);
  Delete(pointer_map, message.in().plugin());
  // TODO(phl): should we do something with out() here?
}

void DirectlyInsertCelestial::Fill(In const& in,
                                   not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
  if (in.parent_index != nullptr) {
    m->set_parent_index(*in.parent_index);
  }
  m->set_gravitational_parameter(in.gravitational_parameter);
  if (in.axis_right_ascension != nullptr) {
    m->set_axis_right_ascension(in.axis_right_ascension);
  }
  if (in.axis_declination != nullptr) {
    m->set_axis_declination(in.axis_declination);
  }
  if (in.j2 != nullptr) {
    m->set_j2(in.j2);
  }
  if (in.reference_radius != nullptr) {
    m->set_reference_radius(in.reference_radius);
  }
  m->set_x(in.x);
  m->set_y(in.y);
  m->set_z(in.z);
  m->set_vx(in.vx);
  m->set_vy(in.vy);
  m->set_vz(in.vz);
}

void DirectlyInsertCelestial::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  int const parent_index = in.parent_index();
  ksp_plugin::principia__DirectlyInsertCelestial(
      plugin,
      in.celestial_index(),
      in.has_parent_index() ? &parent_index : nullptr,
      in.gravitational_parameter().c_str(),
      in.has_axis_right_ascension() ?
          in.axis_right_ascension().c_str() : nullptr,
      in.has_axis_declination() ? in.axis_declination().c_str() : nullptr,
      in.has_j2() ? in.j2().c_str() : nullptr,
      in.has_reference_radius() ? in.reference_radius().c_str() : nullptr,
      in.x().c_str(), in.y().c_str(), in.z().c_str(),
      in.vx().c_str(), in.vy().c_str(), in.vz().c_str());
}

void InsertCelestial::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
  m->set_gravitational_parameter(in.gravitational_parameter);
  m->set_parent_index(in.parent_index);
  *m->mutable_from_parent() = SerializeQP(in.from_parent);
}

void InsertCelestial::Run(Message const& message,
                          not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__InsertCelestial(plugin,
                                         in.celestial_index(),
                                         in.gravitational_parameter(),
                                         in.parent_index(),
                                         DeserializeQP(in.from_parent()));
}

void InsertSun::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
  m->set_gravitational_parameter(in.gravitational_parameter);
}

void InsertSun::Run(Message const& message,
                    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__InsertSun(plugin,
                                   in.celestial_index(),
                                   in.gravitational_parameter());
}

void UpdateCelestialHierarchy::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
  m->set_parent_index(in.parent_index);
}

void UpdateCelestialHierarchy::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__UpdateCelestialHierarchy(plugin,
                                                  in.celestial_index(),
                                                  in.parent_index());
}

void EndInitialization::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void EndInitialization::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map,
                                             message.in().plugin());
  ksp_plugin::principia__EndInitialization(plugin);
}

void InsertOrKeepVessel::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_parent_index(in.parent_index);
}

void InsertOrKeepVessel::Fill(Return const& result,
                              not_null<Message*> const message) {
  auto* m = message->mutable_return_();
  m->set_insert_or_keep_vessel(result);
}

void InsertOrKeepVessel::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().insert_or_keep_vessel(),
           ksp_plugin::principia__InsertOrKeepVessel(plugin,
                                                     in.vessel_guid().c_str(),
                                                     in.parent_index()));
}

void SetVesselStateOffset::Fill(In const& in,
                                not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  *m->mutable_from_parent() = SerializeQP(in.from_parent);
}

void SetVesselStateOffset::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__SetVesselStateOffset(plugin,
                                              in.vessel_guid().c_str(),
                                              DeserializeQP(in.from_parent()));
}

void AdvanceTime::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_t(in.t);
  m->set_planetarium_rotation(in.planetarium_rotation);
}

void AdvanceTime::Run(Message const& message,
                      not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__AdvanceTime(plugin, in.t(), in.planetarium_rotation());
}

void ForgetAllHistoriesBefore::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_t(in.t);
}

void ForgetAllHistoriesBefore::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__ForgetAllHistoriesBefore(plugin, in.t());
}

void VesselFromParent::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void VesselFromParent::Fill(Return const& result,
                            not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_from_parent() =
      SerializeQP(result);
}

void VesselFromParent::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeQP(message.return_().vessel_from_parent()) ==
            ksp_plugin::principia__VesselFromParent(plugin,
                                                    in.vessel_guid().c_str()));
}

void CelestialFromParent::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
}

void CelestialFromParent::Fill(Return const& result,
                               not_null<Message*> const message) {
  *message->mutable_return_()->mutable_celestial_from_parent() =
      SerializeQP(result);
}

void CelestialFromParent::Run(Message const& message,
                              not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeQP(message.return_().celestial_from_parent()) ==
            ksp_plugin::principia__CelestialFromParent(plugin,
                                                       in.celestial_index()));
}

void NewBodyCentredNonRotatingNavigationFrame::Fill(
    In const& in,
    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_reference_body_index(in.reference_body_index);
}

void NewBodyCentredNonRotatingNavigationFrame::Fill(
    Return const& result,
    not_null<Message*> const message) {
  message->mutable_return_()->
      set_new_body_centred_non_rotating_navigation_frame(
          SerializePointer(result));
}

void NewBodyCentredNonRotatingNavigationFrame::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* navigation_frame =
      ksp_plugin::principia__NewBodyCentredNonRotatingNavigationFrame(
          plugin, in.reference_body_index());
  Insert(pointer_map,
         message.return_().new_body_centred_non_rotating_navigation_frame(),
         navigation_frame);
}

void NewBarycentricRotatingNavigationFrame::Fill(
    In const& in,
    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_primary_index(in.primary_index);
  m->set_secondary_index(in.secondary_index);
}

void NewBarycentricRotatingNavigationFrame::Fill(
    Return const& result,
    not_null<Message*> const message) {
  message->mutable_return_()->set_new_barycentric_rotating_navigation_frame(
      SerializePointer(result));
}

void NewBarycentricRotatingNavigationFrame::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* navigation_frame =
      ksp_plugin::principia__NewBarycentricRotatingNavigationFrame(
          plugin, in.primary_index(), in.secondary_index());
  Insert(pointer_map,
         message.return_().new_barycentric_rotating_navigation_frame(),
         navigation_frame);
}

void NewNavigationFrame::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  *m->mutable_parameters() =
      SerializeNavigationFrameParameters(in.parameters);
}

void NewNavigationFrame::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_new_navigation_frame(
      SerializePointer(result));
}

void NewNavigationFrame::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* navigation_frame =
    ksp_plugin::principia__NewNavigationFrame(
        plugin, DeserializeNavigationFrameParameters(in.parameters()));
  Insert(pointer_map,
         message.return_().new_navigation_frame(),
         navigation_frame);
}

void GetNavigationFrameParameters::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_navigation_frame(SerializePointer(in.navigation_frame));
}

void GetNavigationFrameParameters::Fill(Return const& result,
                                        not_null<Message*> const message) {
  *message->mutable_return_()->mutable_get_navigation_frame_parameters() =
      SerializeNavigationFrameParameters(result);
}

void GetNavigationFrameParameters::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* navigation_frame =
      DeserializePointer<NavigationFrame const*>(*pointer_map,
                                                 in.navigation_frame());
  ksp_plugin::principia__GetNavigationFrameParameters(navigation_frame);
  // TODO(phl): Should we check return_() here?
}

void SetPlottingFrame::Fill(In const& in,
                            not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_navigation_frame(SerializePointer(*in.navigation_frame));
}

void SetPlottingFrame::Fill(Out const& out,
                            not_null<Message*> const message) {
  message->mutable_out()->set_navigation_frame(
      SerializePointer(*out.navigation_frame));
}

void SetPlottingFrame::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* navigation_frame = DeserializePointer<NavigationFrame*>(
                               *pointer_map, in.navigation_frame());
  ksp_plugin::principia__SetPlottingFrame(plugin, &navigation_frame);
  // While the |navigation_frame| is not deleted by this function, it is handed
  // over to the plugin so we should not see it cross the interface again.
  Delete(pointer_map, in.navigation_frame());
  // TODO(phl): should we do something with out() here?
}

void GetPlottingFrame::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void GetPlottingFrame::Fill(Return const& result,
                            not_null<Message*> const message) {
  message->mutable_return_()->set_get_plotting_frame(SerializePointer(result));
}

void GetPlottingFrame::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* navigation_frame = ksp_plugin::principia__GetPlottingFrame(plugin);
  Insert(pointer_map,
         message.return_().get_plotting_frame(),
         navigation_frame);

}

void UpdatePrediction::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void UpdatePrediction::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__UpdatePrediction(plugin, in.vessel_guid().c_str());
}

void RenderedVesselTrajectory::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedVesselTrajectory::Fill(Return const& result,
                                    not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_vessel_trajectory(
      SerializePointer(result));
}

void RenderedVesselTrajectory::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* line_and_iterator =
    ksp_plugin::principia__RenderedVesselTrajectory(
        plugin,
        in.vessel_guid().c_str(),
        DeserializeXYZ(in.sun_world_position()));
  Insert(pointer_map,
         message.return_().rendered_vessel_trajectory(),
         line_and_iterator);
}

void HasPrediction::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void HasPrediction::Fill(Return const& result,
                         not_null<Message*> const message) {
  message->mutable_return_()->set_has_prediction(result);
}

void HasPrediction::Run(Message const& message,
                        not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().has_prediction(),
           ksp_plugin::principia__HasPrediction(plugin,
                                                in.vessel_guid().c_str()));
}

void RenderedPrediction::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedPrediction::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_prediction(
      SerializePointer(result));
}

void RenderedPrediction::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* line_and_iterator = ksp_plugin::principia__RenderedPrediction(
                                plugin,
                                in.vessel_guid().c_str(),
                                DeserializeXYZ(in.sun_world_position()));
  Insert(pointer_map,
         message.return_().rendered_prediction(),
         line_and_iterator);
}

void FlightPlanSize::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void FlightPlanSize::Fill(Return const& result,
                          not_null<Message*> const message) {
  message->mutable_return_()->set_flight_plan_size(result);
}

void FlightPlanSize::Run(Message const& message,
                         not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().flight_plan_size(),
           ksp_plugin::principia__FlightPlanSize(plugin,
                                                 in.vessel_guid().c_str()));
}

void RenderedFlightPlan::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_plan_phase(in.plan_phase);
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedFlightPlan::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_flight_plan(
      SerializePointer(result));
}

void RenderedFlightPlan::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* line_and_iterator = ksp_plugin::principia__RenderedFlightPlan(
                                plugin,
                                in.vessel_guid().c_str(),
                                in.plan_phase(),
                                DeserializeXYZ(in.sun_world_position()));
  Insert(pointer_map,
         message.return_().rendered_flight_plan(),
         line_and_iterator);
}

void SetPredictionLength::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_t(in.t);
}

void SetPredictionLength::Run(Message const& message,
                              not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__SetPredictionLength(plugin, in.t());
}

void SetPredictionLengthTolerance::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_l(in.l);
}

void SetPredictionLengthTolerance::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__SetPredictionLengthTolerance(plugin, in.l());
}

void SetPredictionSpeedTolerance::Fill(In const& in,
                                       not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_v(in.v);
}

void SetPredictionSpeedTolerance::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  ksp_plugin::principia__SetPredictionSpeedTolerance(plugin, in.v());
}

void HasVessel::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void HasVessel::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_has_vessel(result);
}

void HasVessel::Run(Message const& message,
                    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().has_vessel(),
           ksp_plugin::principia__HasVessel(plugin, in.vessel_guid().c_str()));
}

void NumberOfSegments::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(in.line_and_iterator));
}

void NumberOfSegments::Fill(Return const& result,
                            not_null<Message*> const message) {
  message->mutable_return_()->set_number_of_segments(result);
}

void NumberOfSegments::Run(Message const& message,
                           not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  CHECK_EQ(message.return_().number_of_segments(),
           ksp_plugin::principia__NumberOfSegments(
               DeserializePointer<LineAndIterator const*>(
                   *pointer_map, in.line_and_iterator())));
}

void FetchAndIncrement::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(in.line_and_iterator));
}

void FetchAndIncrement::Fill(Return const& result,
                             not_null<Message*> const message) {
  *message->mutable_return_()->mutable_fetch_and_increment() =
      SerializeXYZSegment(result);
}

void FetchAndIncrement::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  CHECK(DeserializeXYZSegment(message.return_().fetch_and_increment()) ==
            ksp_plugin::principia__FetchAndIncrement(
                DeserializePointer<LineAndIterator*>(
                    *pointer_map, in.line_and_iterator())));
}

void AtEnd::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(in.line_and_iterator));
}

void AtEnd::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_at_end(result);
}

void AtEnd::Run(Message const& message,
                not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  CHECK_EQ(message.return_().at_end(),
           ksp_plugin::principia__AtEnd(DeserializePointer<LineAndIterator*>(
               *pointer_map, in.line_and_iterator())));
}

void DeleteLineAndIterator::Fill(In const& in,
                                 not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(*in.line_and_iterator));
}

void DeleteLineAndIterator::Fill(Out const& out,
                                 not_null<Message*> const message) {
  message->mutable_out()->set_line_and_iterator(
      SerializePointer(*out.line_and_iterator));
}

void DeleteLineAndIterator::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto* line_and_iterator = DeserializePointer<LineAndIterator*>(
                                *pointer_map, message.in().line_and_iterator());
  ksp_plugin::principia__DeleteLineAndIterator(&line_and_iterator);
  Delete(pointer_map, message.in().line_and_iterator());
  // TODO(phl): should we do something with out() here?
}

void AddVesselToNextPhysicsBubble::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  for (KSPPart const* part = in.parts; part < in.parts + in.count; ++part) {
    *m->add_parts() = SerializeKSPPart(*part);
  }
}

void AddVesselToNextPhysicsBubble::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  std::vector<KSPPart> deserialized_parts;
  deserialized_parts.reserve(in.parts_size());
  for (auto const& part : in.parts()) {
    deserialized_parts.push_back(DeserializeKSPPart(part));
  }
  ksp_plugin::principia__AddVesselToNextPhysicsBubble(
      plugin,
      in.vessel_guid().c_str(),
      &deserialized_parts[0],
      deserialized_parts.size());
}

void PhysicsBubbleIsEmpty::Fill(In const& in,
                                not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void PhysicsBubbleIsEmpty::Fill(Return const& result,
                                not_null<Message*> const message) {
  message->mutable_return_()->set_physics_bubble_is_empty(result);
}

void PhysicsBubbleIsEmpty::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().physics_bubble_is_empty(),
           ksp_plugin::principia__PhysicsBubbleIsEmpty(plugin));
}

void BubbleDisplacementCorrection::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  *m->mutable_sun_position() = SerializeXYZ(in.sun_position);
}

void BubbleDisplacementCorrection::Fill(Return const& result,
                                        not_null<Message*> const message) {
  *message->mutable_return_()->mutable_bubble_displacement_correction() =
      SerializeXYZ(result);
}

void BubbleDisplacementCorrection::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeXYZ(message.return_().bubble_displacement_correction()) ==
            ksp_plugin::principia__BubbleDisplacementCorrection(
                plugin, DeserializeXYZ(in.sun_position())));
}

void BubbleVelocityCorrection::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_reference_body_index(in.reference_body_index);
}

void BubbleVelocityCorrection::Fill(Return const& result,
                                    not_null<Message*> const message) {
  *message->mutable_return_()->mutable_bubble_velocity_correction() =
      SerializeXYZ(result);
}

void BubbleVelocityCorrection::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeXYZ(message.return_().bubble_velocity_correction()) ==
            ksp_plugin::principia__BubbleVelocityCorrection(
                plugin, in.reference_body_index()));
}

void NavballOrientation::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
  *m->mutable_ship_world_position() = SerializeXYZ(in.ship_world_position);
}

void NavballOrientation::Fill(Return const& result,
                              not_null<Message*> const message) {
  *message->mutable_return_()->mutable_navball_orientation() =
      SerializeWXYZ(result);
}

void NavballOrientation::Run(Message const& message,
                             not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeWXYZ(message.return_().navball_orientation()) ==
            ksp_plugin::principia__NavballOrientation(
                plugin,
                DeserializeXYZ(in.sun_world_position()),
                DeserializeXYZ(in.ship_world_position())));
}

void VesselTangent::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void VesselTangent::Fill(Return const& result,
                         not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_tangent() =
      SerializeXYZ(result);
}

void VesselTangent::Run(Message const& message,
                        not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeXYZ(message.return_().vessel_tangent()) ==
            ksp_plugin::principia__VesselTangent(
                plugin,
                in.vessel_guid().c_str()));
}

void VesselNormal::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void VesselNormal::Fill(Return const& result,
                        not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_normal() =
      SerializeXYZ(result);
}

void VesselNormal::Run(Message const& message,
                       not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeXYZ(message.return_().vessel_normal()) ==
            ksp_plugin::principia__VesselNormal(
                plugin,
                in.vessel_guid().c_str()));
}

void VesselBinormal::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void VesselBinormal::Fill(Return const& result,
                          not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_binormal() =
      SerializeXYZ(result);
}

void VesselBinormal::Run(Message const& message,
                         not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK(DeserializeXYZ(message.return_().vessel_binormal()) ==
            ksp_plugin::principia__VesselBinormal(
                plugin,
                in.vessel_guid().c_str()));
}

void CurrentTime::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void CurrentTime::Fill(Return const& result,
                       not_null<Message*> const message) {
  message->mutable_return_()->set_current_time(result);
}

void CurrentTime::Run(Message const& message,
                      not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  CHECK_EQ(message.return_().current_time(),
           ksp_plugin::principia__CurrentTime(plugin));
}

void SerializePlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_serializer(SerializePointer(*in.serializer));
}

void SerializePlugin::Fill(Out const& out, not_null<Message*> const message) {
  auto* m = message->mutable_out();
  m->set_serializer(SerializePointer(*out.serializer));
}

void SerializePlugin::Fill(Return const& result,
                           not_null<Message*> const message) {
  if (result != nullptr) {
    message->mutable_return_()->set_serialize_plugin(SerializePointer(result));
  }
}

void SerializePlugin::Run(Message const& message,
                          not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* serializer = DeserializePointer<PullSerializer*>(
                         *pointer_map, in.serializer());
  char const* serialize_plugin =
      ksp_plugin::principia__SerializePlugin(plugin, &serializer);
  if (serialize_plugin == nullptr) {
    Delete(pointer_map, in.serializer());
  } else {
    Insert(pointer_map, message.out().serializer(), serializer);
    Insert(pointer_map, message.return_().serialize_plugin(), serialize_plugin);
  }
}

void DeletePluginSerialization::Fill(In const& in,
                                     not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_serialization(SerializePointer(*in.serialization));
}

void DeletePluginSerialization::Fill(Out const& out,
                                     not_null<Message*> const message) {
  auto* m = message->mutable_out();
  m->set_serialization(SerializePointer(*out.serialization));
}

void DeletePluginSerialization::Run(
    Message const& message,
    not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* serialization = DeserializePointer<char const*>(
                            *pointer_map, in.serialization());
  ksp_plugin::principia__DeletePluginSerialization(&serialization);
  Delete(pointer_map, in.serialization());
}

void DeserializePlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_serialization(std::string(in.serialization, in.serialization_size));
  m->set_deserializer(SerializePointer(*in.deserializer));
  m->set_plugin(SerializePointer(*in.plugin));
}

void DeserializePlugin::Fill(Out const& out, not_null<Message*> const message) {
  auto* m = message->mutable_out();
  m->set_deserializer(SerializePointer(*out.deserializer));
  m->set_plugin(SerializePointer(*out.plugin));
}

void DeserializePlugin::Run(Message const& message,
                            not_null<Player::PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin const*>(*pointer_map, in.plugin());
  auto* deserializer = DeserializePointer<PushDeserializer*>(
                           *pointer_map, in.deserializer());
  ksp_plugin::principia__DeserializePlugin(in.serialization().c_str(),
                                           in.serialization().size(),
                                           &deserializer,
                                           &plugin);
  if (in.serialization().size() == 0) {
    Delete(pointer_map, in.deserializer());
  } else {
    Insert(pointer_map, message.out().deserializer(), deserializer);
  }
  Insert(pointer_map, message.out().plugin(), plugin);
}

void SayHello::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_say_hello(result);
}

void SayHello::Run(Message const& message,
                   not_null<Player::PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().say_hello(),
           ksp_plugin::principia__SayHello());
}

}  // namespace journal
}  // namespace principia
