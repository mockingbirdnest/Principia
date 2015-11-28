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

int const kBufferSize = 100;

template<typename T,
         typename = typename std::enable_if<std::is_pointer<T>::value>::type>
T Find(PointerMap const& pointer_map, std::uint64_t const address) {
    return reinterpret_cast<T>(FindOrDie(pointer_map, address));
}

template<typename T>
void Insert(not_null<PointerMap*> const pointer_map,
            std::uint64_t const address,
            T* const pointer) {
  auto inserted = pointer_map->emplace(address, pointer);
  CHECK(inserted.second) << address;
}

XYZ DeserializeXYZ(serialization::XYZ const& xyz) {
  return {xyz.x(), xyz.y(), xyz.z()};
}

QP DeserializeQP(serialization::QP const& qp) {
  return {DeserializeXYZ(qp.q()), DeserializeXYZ(qp.p())};
}

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

serialization::XYZ SerializeXYZ(XYZ const& xyz) {
  serialization::XYZ m;
  m.set_x(xyz.x);
  m.set_y(xyz.y);
  m.set_z(xyz.z);
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
  message->mutable_return_()->set_max_severity(result);
}

void GetBufferedLogging::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().max_severity(), principia__GetBufferedLogging());
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
  message->mutable_return_()->set_seconds(result);
}

void GetBufferDuration::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().seconds(), principia__GetBufferDuration());
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
  message->mutable_return_()->set_min_severity(result);
}

void GetSuppressedLogging::Run(Message const& message,
                               not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().min_severity(), principia__GetSuppressedLogging());
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
  message->mutable_return_()->set_level(result);
}

void GetVerboseLogging::Run(Message const& message,
                            not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().level(), principia__GetVerboseLogging());
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
  message->mutable_return_()->set_min_severity(result);
}

void GetStderrLogging::Run(Message const& message,
                           not_null<PointerMap*> const pointer_map) {
  CHECK_EQ(message.return_().min_severity(), principia__GetStderrLogging());
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
  message->mutable_return_()->set_plugin(SerializePointer(result));
}

void NewPlugin::Run(Message const& message,
                    not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = principia__NewPlugin(in.initial_time(),
                                      in.planetarium_rotation_in_degrees());
  Insert(pointer_map, message.return_().plugin(), plugin);
}

void DeletePlugin::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_plugin(SerializePointer(in.plugin));
}

void DeletePlugin::Fill(Out const& out, not_null<Message*> const message) {
  message->mutable_out()->set_plugin(SerializePointer(*out.plugin));
}

void DeletePlugin::Run(Message const& message,
                       not_null<PointerMap*> const pointer_map) {
  auto* plugin = Find<Plugin const*>(*pointer_map, message.in().plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, message.in().plugin());
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
  m->set_inserted(result);
}

void InsertOrKeepVessel::Run(Message const& message,
                             not_null<PointerMap*> const pointer_map) {
  auto const& in = message.in();
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
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
  auto* plugin = Find<Plugin*>(*pointer_map, in.plugin());
  principia__ForgetAllHistoriesBefore(plugin, in.t());
}

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
