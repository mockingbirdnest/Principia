#include "ksp_plugin/journal.hpp"

#include <fstream>
#include <list>
#include <string>
#include <type_traits>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "base/map_util.hpp"
#include "glog/logging.h"

namespace principia {
namespace ksp_plugin {

using base::Bytes;
using base::FindOrDie;
using base::GetLine;
using base::HexadecimalDecode;
using base::HexadecimalEncode;
using base::UniqueBytes;

namespace {

template<typename T>
void Insert(not_null<PointerMap*> const pointer_map,
            std::uint64_t const address,
            T* const pointer) {
  auto inserted = pointer_map->emplace(address, pointer);
  CHECK(inserted.second) << address;
}

template<typename T,
         typename = typename std::enable_if<std::is_pointer<T>::value>::type>
T DeserializePointer(PointerMap const& pointer_map,
                     std::uint64_t const address) {
  return reinterpret_cast<T>(FindOrDie(pointer_map, address));
}

XYZ DeserializeXYZ(serialization::XYZ const& xyz) {
  return {xyz.x(), xyz.y(), xyz.z()};
}

QP DeserializeQP(serialization::QP const& qp) {
  return {DeserializeXYZ(qp.q()), DeserializeXYZ(qp.p())};
}

serialization::KSPPart SerializePart(KSPPart const& ksp_part) {
  serialization::KSPPart m;
  *m.mutable_world_position() = SerializeXYZ(ksp_part.world_position);
  *m.mutable_world_velocity() = SerializeXYZ(ksp_part.world_velocity);
  m.set_mass(ksp_part.mass);
  *m.mutable_gravitational_acceleration_to_be_applied_by_ksp() =
      SerializeXYZ(ksp_part.gravitational_acceleration_to_be_applied_by_ksp);
  m.set_id(ksp_part.id);
}

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

serialization::WXYZ SerializeWXYZ(WXYZ const& wxyz) {
  serialization::WXYZ m;
  m.set_x(wxyz.w);
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

}  // namespace

void SetBufferedLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_max_severity(in.max_severity);
}

void SetBufferedLogging::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  principia__SetBufferedLogging(message.in().max_severity());
}

void GetBufferedLogging::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_get_buffered_logging(result);
}

void GetBufferedLogging::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_buffered_logging(),
           principia__GetBufferedLogging());
}

void SetBufferDuration::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_seconds(in.seconds);
}

void SetBufferDuration::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  principia__SetBufferDuration(message.in().seconds());
}

void GetBufferDuration::Fill(Return const& result,
                             not_null<Message*> const message) {
  message->mutable_return_()->set_get_buffer_duration(result);
}

void GetBufferDuration::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_buffer_duration(),
           principia__GetBufferDuration());
}

void SetSuppressedLogging::Fill(In const& in,
                                not_null<Message*> const message) {
  message->mutable_in()->set_min_severity(in.min_severity);
}

void SetSuppressedLogging::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {
  principia__SetSuppressedLogging(message.in().min_severity());
}

void GetSuppressedLogging::Fill(Return const& result,
                                not_null<Message*> const message) {
  message->mutable_return_()->set_get_suppressed_logging(result);
}

void GetSuppressedLogging::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_suppressed_logging(),
           principia__GetSuppressedLogging());
}

void SetVerboseLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_level(in.level);
}

void SetVerboseLogging::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  principia__SetVerboseLogging(message.in().level());
}

void GetVerboseLogging::Fill(Return const& result,
                             not_null<Message*> const message) {
  message->mutable_return_()->set_get_verbose_logging(result);
}

void GetVerboseLogging::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_verbose_logging(),
           principia__GetVerboseLogging());
}

void SetStderrLogging::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_min_severity(in.min_severity);
}

void SetStderrLogging::Run(Message const& message,
                           not_null<PointerMap*> const pointer_map) {
  principia__SetStderrLogging(message.in().min_severity());
}

void GetStderrLogging::Fill(Return const& result,
                            not_null<Message*> const message) {
  message->mutable_return_()->set_get_stderr_logging(result);
}

void GetStderrLogging::Run(Message const& message,
                           not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().get_stderr_logging(),
           principia__GetStderrLogging());
}

void LogInfo::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogInfo::Run(Message const& message,
                  not_null<PointerMap*> const pointer_map) {
  principia__LogInfo(message.in().message().c_str());
}

void LogWarning::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogWarning::Run(Message const& message,
                     not_null<PointerMap*> const pointer_map) {
  principia__LogWarning(message.in().message().c_str());
}

void LogError::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogError::Run(Message const& message,
                   not_null<PointerMap*> const pointer_map) {
  principia__LogError(message.in().message().c_str());
}

void LogFatal::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_message(in.message);
}

void LogFatal::Run(Message const& message,
                   not_null<PointerMap*> const pointer_map) {
  principia__LogFatal(message.in().message().c_str());
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
                    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = principia__NewPlugin(in.initial_time(),
                                      in.planetarium_rotation_in_degrees());
  Insert(pointer_map, message.return_().new_plugin(), plugin);
}

void DeletePlugin::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_plugin(SerializePointer(in.plugin));
}

void DeletePlugin::Fill(Out const& out, not_null<Message*> const message) {
  message->mutable_out()->set_plugin(SerializePointer(*out.plugin));
}

void DeletePlugin::Run(Message const& message,
                       not_null<PointerMap*> const pointer_map) {
  auto* plugin = DeserializePointer<Plugin const*>(*pointer_map,
                                                   message.in().plugin());
  principia__DeletePlugin(&plugin);
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

void DirectlyInsertCelestial::Run(Message const& message,
                                  not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  int const parent_index = in.parent_index();
  principia__DirectlyInsertCelestial(
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
                          not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__InsertCelestial(plugin,
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
                    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__InsertSun(plugin,
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

void UpdateCelestialHierarchy::Run(Message const& message,
                                   not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__UpdateCelestialHierarchy(plugin,
                                      in.celestial_index(),
                                      in.parent_index());
}

void EndInitialization::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void EndInitialization::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map,
                                             message.in().plugin());
  principia__EndInitialization(plugin);
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
                             not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__InsertOrKeepVessel(plugin,
                                in.vessel_guid().c_str(),
                                in.parent_index());
  // TODO(phl): should we do something with out() here?
}

void SetVesselStateOffset::Fill(In const& in,
                                not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  *m->mutable_from_parent() = SerializeQP(in.from_parent);
}

void SetVesselStateOffset::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__SetVesselStateOffset(plugin,
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
                      not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__AdvanceTime(plugin, in.t(), in.planetarium_rotation());
}

void ForgetAllHistoriesBefore::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_t(in.t);
}

void ForgetAllHistoriesBefore::Run(Message const& message,
                                   not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__ForgetAllHistoriesBefore(plugin, in.t());
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
                           not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__VesselFromParent(plugin, in.vessel_guid().c_str());
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
                              not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__CelestialFromParent(plugin, in.celestial_index());
  // TODO(phl): Check all the return values everywhere.
}

void NewBodyCentredNonRotatingRenderingFrame::Fill(
    In const& in,
    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_reference_body_index(in.reference_body_index);
}

void NewBodyCentredNonRotatingRenderingFrame::Fill(
    Return const& result,
    not_null<Message*> const message) {
  message->mutable_return_()->set_new_body_centred_non_rotating_rendering_frame(
      SerializePointer(result));
}

void NewBodyCentredNonRotatingRenderingFrame::Run(
    Message const& message,
    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* rendering_frame = principia__NewBodyCentredNonRotatingRenderingFrame(
                              plugin, in.reference_body_index());
  Insert(pointer_map,
         message.return_().new_body_centred_non_rotating_rendering_frame(),
         rendering_frame);
}

void NewBarycentricRotatingRenderingFrame::Fill(
    In const& in,
    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_primary_index(in.primary_index);
  m->set_secondary_index(in.secondary_index);
}

void NewBarycentricRotatingRenderingFrame::Fill(
    Return const& result,
    not_null<Message*> const message) {
  message->mutable_return_()->set_new_barycentric_rotating_rendering_frame(
      SerializePointer(result));
}

void NewBarycentricRotatingRenderingFrame::Run(
    Message const& message,
    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* rendering_frame = principia__NewBarycentricRotatingRenderingFrame(
                              plugin, in.primary_index(), in.secondary_index());
  Insert(pointer_map,
         message.return_().new_barycentric_rotating_rendering_frame(),
         rendering_frame);
}

void DeleteRenderingFrame::Fill(In const& in,
                                not_null<Message*> const message) {
  message->mutable_in()->set_rendering_frame(
      SerializePointer(in.rendering_frame));
}

void DeleteRenderingFrame::Fill(Out const& out,
                                not_null<Message*> const message) {
  message->mutable_out()->set_rendering_frame(
      SerializePointer(*out.rendering_frame));
}

void DeleteRenderingFrame::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {
  auto* rendering_frame = DeserializePointer<RenderingFrame*>(
                              *pointer_map, message.in().rendering_frame());
  principia__DeleteRenderingFrame(&rendering_frame);
  // TODO(phl): should we do something with out() here?
}

void UpdatePrediction::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
}

void UpdatePrediction::Run(Message const& message,
                           not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__UpdatePrediction(plugin, in.vessel_guid().c_str());
}

void RenderedVesselTrajectory::Fill(In const& in,
                                    not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedVesselTrajectory::Fill(Return const& result,
                                    not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_vessel_trajectory(
      SerializePointer(result));
}

void RenderedVesselTrajectory::Run(Message const& message,
                                   not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* rendering_frame = DeserializePointer<RenderingFrame*>(
                              *pointer_map, in.rendering_frame());
  auto* line_and_iterator = principia__RenderedVesselTrajectory(
                                plugin,
                                in.vessel_guid().c_str(),
                                rendering_frame,
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
                        not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__HasPrediction(plugin, in.vessel_guid().c_str());
}

void RenderedPrediction::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedPrediction::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_prediction(
      SerializePointer(result));
}

void RenderedPrediction::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* rendering_frame = DeserializePointer<RenderingFrame*>(
                              *pointer_map, in.rendering_frame());
  auto* line_and_iterator = principia__RenderedPrediction(
                                plugin,
                                in.vessel_guid().c_str(),
                                rendering_frame,
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
                         not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__FlightPlanSize(plugin, in.vessel_guid().c_str());
}

void RenderedFlightPlan::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_plan_phase(in.plan_phase);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
  *m->mutable_sun_world_position() = SerializeXYZ(in.sun_world_position);
}

void RenderedFlightPlan::Fill(Return const& result,
                              not_null<Message*> const message) {
  message->mutable_return_()->set_rendered_flight_plan(
      SerializePointer(result));
}

void RenderedFlightPlan::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  auto* line_and_iterator = principia__RenderedFlightPlan(
                                plugin,
                                in.vessel_guid().c_str(),
                                in.plan_phase(),
                                DeserializePointer<RenderingFrame*>(
                                    *pointer_map, in.rendering_frame()),
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
                              not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__SetPredictionLength(plugin, in.t());
}

void SetPredictionLengthTolerance::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_l(in.l);
}

void SetPredictionLengthTolerance::Run(
    Message const& message,
    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__SetPredictionLengthTolerance(plugin, in.l());
}

void SetPredictionSpeedTolerance::Fill(In const& in,
                                       not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_v(in.v);
}

void SetPredictionSpeedTolerance::Run(Message const& message,
                                      not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__SetPredictionSpeedTolerance(plugin, in.v());
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
                    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = DeserializePointer<Plugin*>(*pointer_map, in.plugin());
  principia__HasVessel(plugin, in.vessel_guid().c_str());
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
                           not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  principia__NumberOfSegments(DeserializePointer<LineAndIterator const*>(
                                  *pointer_map, in.line_and_iterator()));
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
                            not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  principia__FetchAndIncrement(DeserializePointer<LineAndIterator*>(
                                   *pointer_map, in.line_and_iterator()));
}

void AtEnd::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(in.line_and_iterator));
}

void AtEnd::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_at_end(result);
}

void AtEnd::Run(Message const& message,
                not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  principia__AtEnd(DeserializePointer<LineAndIterator*>(
                       *pointer_map, in.line_and_iterator()));
}

void DeleteLineAndIterator::Fill(In const& in,
                                 not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_line_and_iterator(SerializePointer(in.line_and_iterator));
}

void DeleteLineAndIterator::Fill(Out const& out,
                                 not_null<Message*> const message) {
  message->mutable_out()->set_line_and_iterator(
      SerializePointer(*out.line_and_iterator));
}

void DeleteLineAndIterator::Run(Message const& message,
                                not_null<PointerMap*> const pointer_map) {
  auto* line_and_iterator = DeserializePointer<LineAndIterator*>(
                                *pointer_map, message.in().line_and_iterator());
  principia__DeleteLineAndIterator(&line_and_iterator);
  // TODO(phl): should we do something with out() here?
}

void AddVesselToNextPhysicsBubble::Fill(In const& in,
                                        not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  for (KSPPart const* part = in.parts; part < in.parts + in.count; ++part) {
    *m->add_parts() = SerializePart(*part);
  }
}

void AddVesselToNextPhysicsBubble::Run(
    Message const& message,
    not_null<PointerMap*> const pointer_map) {}

void PhysicsBubbleIsEmpty::Fill(In const& in,
                                not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void PhysicsBubbleIsEmpty::Fill(Return const& result,
                                not_null<Message*> const message) {
  message->mutable_return_()->set_physics_buble_is_empty(result);
}

void PhysicsBubbleIsEmpty::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {}

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

void BubbleDisplacementCorrection::Run(Message const& message,
                                       not_null<PointerMap*> const pointer_map) {}

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

void BubbleVelocityCorrection::Run(Message const& message,
                                   not_null<PointerMap*> const pointer_map) {}

void NavballOrientation::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void NavballOrientation::Fill(Return const& result,
                              not_null<Message*> const message) {
  *message->mutable_return_()->mutable_navball_orientation() =
      SerializeWXYZ(result);
}

void NavballOrientation::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {}

void VesselTangent::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
}

void VesselTangent::Fill(Return const& result,
                         not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_tangent() =
      SerializeXYZ(result);
}

void VesselTangent::Run(Message const& message,
                        not_null<PointerMap*> const pointer_map) {}

void VesselNormal::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
}

void VesselNormal::Fill(Return const& result,
                        not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_normal() =
      SerializeXYZ(result);
}

void VesselNormal::Run(Message const& message,
                       not_null<PointerMap*> const pointer_map) {}

void VesselBinormal::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_vessel_guid(in.vessel_guid);
  m->set_rendering_frame(SerializePointer(in.rendering_frame));
}

void VesselBinormal::Fill(Return const& result,
                          not_null<Message*> const message) {
  *message->mutable_return_()->mutable_vessel_binormal() =
      SerializeXYZ(result);
}

void VesselBinormal::Run(Message const& message,
                         not_null<PointerMap*> const pointer_map) {}

void CurrentTime::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
}

void CurrentTime::Fill(Return const& result,
                       not_null<Message*> const message) {
  message->mutable_return_()->set_current_time(result);
}

void CurrentTime::Run(Message const& message,
                      not_null<PointerMap*> const pointer_map) {}

void SerializePlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_serializer(SerializePointer(in.serializer));
}

void SerializePlugin::Fill(Out const& out, not_null<Message*> const message) {}

void SerializePlugin::Fill(Return const& result,
                           not_null<Message*> const message) {
  message->mutable_return_()->set_serialize_plugin(result);
}

void SerializePlugin::Run(Message const& message,
                          not_null<PointerMap*> const pointer_map) {}

void DeletePluginSerialization::Fill(In const& in,
                                     not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_serialization(SerializePointer(in.serialization));
}

void DeletePluginSerialization::Fill(Out const& out,
                                     not_null<Message*> const message) {}

void DeletePluginSerialization::Run(Message const& message,
                                    not_null<PointerMap*> const pointer_map) {}

void DeserializePlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* m = message->mutable_in();
  m->set_serialization(std::string(in.serialization, in.serialization_size));
  m->set_deserializer(SerializePointer(in.deserializer));
  m->set_plugin(SerializePointer(in.plugin));
}

void DeserializePlugin::Fill(Out const& out, not_null<Message*> const message) {}

void DeserializePlugin::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {}

void SayHello::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_say_hello(result);
}

void SayHello::Run(Message const& message,
                   not_null<PointerMap*> const pointer_map) {}

Journal::Journal(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::out) {}

Journal::~Journal() {
  stream_.close();
}

void Journal::Write(serialization::Method const& method) {
  UniqueBytes bytes(method.ByteSize());
  method.SerializeToArray(bytes.data.get(), static_cast<int>(bytes.size));

  std::int64_t const hexadecimal_size = (bytes.size << 1) + 2;
  UniqueBytes hexadecimal(hexadecimal_size);
  HexadecimalEncode({bytes.data.get(), bytes.size}, hexadecimal.get());
  hexadecimal.data.get()[hexadecimal_size - 2] = '\n';
  hexadecimal.data.get()[hexadecimal_size - 1] = '\0';
  stream_ << hexadecimal.data.get();
  stream_.flush();
}

void Journal::Activate(base::not_null<Journal*> const journal) {
  CHECK(active_ == nullptr);
  active_ = journal;
}

void Journal::Deactivate() {
  CHECK(active_ != nullptr);
  delete active_;
  active_ = nullptr;
}

Player::Player(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::in) {}

bool Player::Play() {
  std::unique_ptr<serialization::Method> method = Read();
  if (method == nullptr) {
    return false;
  }

  bool ran = false;
  ran |= RunIfAppropriate<DeletePlugin>(*method);
  ran |= RunIfAppropriate<NewPlugin>(*method);
  CHECK(ran);

  return true;
}

std::unique_ptr<serialization::Method> Player::Read() {
  std::string const line = GetLine(&stream_);
  if (line.empty()) {
    return nullptr;
  }

  uint8_t const* const hexadecimal =
      reinterpret_cast<uint8_t const*>(line.c_str());
  int const hexadecimal_size = strlen(line.c_str());
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(),
                               static_cast<int>(bytes.size)));

  return method;
}

Journal* Journal::active_ = nullptr;

}  // namespace ksp_plugin
}  // namespace principia
