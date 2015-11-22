#include "ksp_plugin/journal.hpp"

#include <fstream>
#include <list>
#include <string>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "glog/logging.h"

namespace principia {
namespace ksp_plugin {

using base::Bytes;
using base::HexadecimalEncode;
using base::UniqueBytes;

namespace {

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

void NewPlugin::Fill(In const& in, not_null<Message*> const message) {
  auto* mutable_in = message->mutable_in();
  mutable_in->set_initial_time(in.initial_time);
  mutable_in->set_planetarium_rotation_in_degrees(
      in.planetarium_rotation_in_degrees);
}

void NewPlugin::Fill(Return const& result, not_null<Message*> const message) {
  message->mutable_return_()->set_plugin(SerializePointer(result));
}

void DeletePlugin::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_plugin(SerializePointer(in.plugin));
}

void DeletePlugin::Fill(Out const& out, not_null<Message*> const message) {
  message->mutable_out()->set_plugin(SerializePointer(*out.plugin));
}

void DirectlyInsertCelestial::Fill(In const& in,
                                   not_null<Message*> const message) {
  Message::In* m = message->mutable_in();
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

void InsertCelestial::Fill(In const& in, not_null<Message*> const message) {
  Message::In* m = message->mutable_in();
  m->set_plugin(SerializePointer(in.plugin));
  m->set_celestial_index(in.celestial_index);
  m->set_gravitational_parameter(in.gravitational_parameter);
  m->set_parent_index(in.parent_index);
  *m->mutable_from_parent() = SerializeQP(in.from_parent);
}

Journal::Journal(std::string const& filename)
    : stream_(filename.c_str(), std::ios::out) {}

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

Journal* Journal::active_ = nullptr;

}  // namespace ksp_plugin
}  // namespace principia
