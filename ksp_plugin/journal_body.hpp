#pragma once

#include "glog/logging.h"
#include "ksp_plugin/journal.hpp"

namespace principia {
namespace ksp_plugin {

template<typename Profile>
Journal::Method<Profile>::Method(typename Profile::In const& in,
                                 typename Profile::Out const& out)
    : message_(std::make_unique<Profile::Message>()),
      out_(out) {
  Profile::Fill(in, message_.get());
}

template<typename Profile>
Journal::Method<Profile>::Method(typename Profile::In const& in)
    : message_(std::make_unique<typename Profile::Message>()) {
  Profile::Fill(in, message_.get());
}

template<typename Profile>
Journal::Method<Profile>::~Method() {
  CHECK(returned_);
  if (out_) {
    Profile::Fill(*out_, message_.get());
  }
}

template<typename Profile>
typename Profile::Return Journal::Method<Profile>::Return(
    typename Profile::Return const& r3turn) {
  CHECK(!returned_);
  returned_ = true;
  Profile::Fill(r3turn, message_.get());
  return r3turn;
}

template<typename Profile>
void Journal::Method<Profile>::Return() {
  CHECK(!returned_);
  returned_ = true;
}

}  // namespace ksp_plugin
}  // namespace principia
