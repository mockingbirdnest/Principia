#pragma once

#include "glog/logging.h"
#include "ksp_plugin/journal.hpp"

namespace principia {
namespace ksp_plugin {

namespace {

template<typename T>
std::uint64_t SerializePointer(T* t) {
  return reinterpret_cast<std::uint64_t>(t);
}

}  // namespace

void DeletePlugin::Fill(In const& in, not_null<Message*> const message) {
  message->mutable_in()->set_plugin(SerializePointer(in.plugin));
}

void DeletePlugin::Fill(Out const& out, not_null<Message*> const message) {
  message->mutable_out()->set_plugin(SerializePointer(*out.plugin));
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

template<typename Profile>
Journal::Method<Profile>::Method(typename Profile::In const& in)
    : message_(std::make_unique<typename Profile::Message>()) {
  Profile::Fill(in, message_.get());
}

template<typename Profile>
template<typename P, typename>
Journal::Method<Profile>::Method(typename Profile::In const& in,
                                 typename P::Out const& out)
    : message_(std::make_unique<typename Profile::Message>()),
      out_filler_([this, out]() { Profile::Fill(out, message_.get()); }) {
  Profile::Fill(in, message_.get());
}

template<typename Profile>
Journal::Method<Profile>::~Method() {
  CHECK(returned_);
  if (out_filler_ != nullptr) {
    out_filler_();
  }
}

template<typename Profile>
void Journal::Method<Profile>::Return() {
  CHECK(!returned_);
  returned_ = true;
}

template<typename Profile>
template<typename P, typename>
typename P::Return Journal::Method<Profile>::Return(
    typename P::Return const& result) {
  CHECK(!returned_);
  returned_ = true;
  Profile::Fill(result, message_.get());
  return result;
}

}  // namespace ksp_plugin
}  // namespace principia
