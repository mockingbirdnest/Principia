#pragma once

#include "journal/method.hpp"

#include <list>

#include "glog/logging.h"
#include "journal/recorder.hpp"

namespace principia {
namespace journal {

template<typename Profile>
Method<Profile>::Method() {
  if (Recorder::active_recorder_ != nullptr) {
    message_ = std::make_unique<typename Profile::Message>();
  }
}

template<typename Profile>
template<typename P, typename>
Method<Profile>::Method(typename P::In const& in) {
  if (Recorder::active_recorder_ != nullptr) {
    message_ = std::make_unique<typename Profile::Message>();
    Profile::Fill(in, message_.get());
  }
}

template<typename Profile>
template<typename P, typename, typename>
Method<Profile>::Method(typename P::In const& in, typename P::Out const& out) {
  if (Recorder::active_recorder_ != nullptr) {
    message_ = std::make_unique<typename Profile::Message>();
    out_filler_ = [this, out]() { Profile::Fill(out, message_.get()); };
    Profile::Fill(in, message_.get());
  }
}

template<typename Profile>
Method<Profile>::~Method() {
  CHECK(returned_);
  if (Recorder::active_recorder_ != nullptr) {
    if (out_filler_ != nullptr) {
      out_filler_();
    }
    serialization::Method method;
    method.SetAllocatedExtension(
        Profile::Message::extension, message_.release());
    Recorder::active_recorder_->Write(method);
  }
}

template<typename Profile>
void Method<Profile>::Return() {
  CHECK(!returned_);
  returned_ = true;
}

template<typename Profile>
template<typename P, typename>
typename P::Return Method<Profile>::Return(
    typename P::Return const& result) {
  CHECK(!returned_);
  returned_ = true;
  if (Recorder::active_recorder_ != nullptr) {
    Profile::Fill(result, message_.get());
  }
  return result;
}

}  // namespace journal
}  // namespace principia
