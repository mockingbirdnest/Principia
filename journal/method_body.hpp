
#pragma once

#include "journal/method.hpp"

#include <list>

#include "journal/recorder.hpp"

namespace principia {
namespace journal {
namespace internal_method {

template<typename Profile>
Method<Profile>::Method() {
  if (Recorder::active_recorder_ != nullptr) {
    serialization::Method method;
    [[maybe_unused]] auto* const message_in =
        method.MutableExtension(Profile::Message::extension);
    Recorder::active_recorder_->WriteAtConstruction(method);
  }
}

template<typename Profile>
template<typename P, typename>
Method<Profile>::Method(typename P::In const& in) {
  if (Recorder::active_recorder_ != nullptr) {
    serialization::Method method;
    auto* const message_in =
        method.MutableExtension(Profile::Message::extension);
    Profile::Fill(in, message_in);
    Recorder::active_recorder_->WriteAtConstruction(method);
  }
}

template<typename Profile>
template<typename P, typename>
Method<Profile>::Method(typename P::Out const& out) {
  if (Recorder::active_recorder_ != nullptr) {
    serialization::Method method;
    [[maybe_unused]] auto* const message_in =
        method.MutableExtension(Profile::Message::extension);
    Recorder::active_recorder_->WriteAtConstruction(method);
    out_filler_ = [out](
        not_null<typename Profile::Message*> const message) {
      Profile::Fill(out, message);
    };
  }
}

template<typename Profile>
template<typename P, typename>
Method<Profile>::Method(typename P::In const& in, typename P::Out const& out) {
  if (Recorder::active_recorder_ != nullptr) {
    serialization::Method method;
    auto* const message_in =
        method.MutableExtension(Profile::Message::extension);
    Profile::Fill(in, message_in);
    Recorder::active_recorder_->WriteAtConstruction(method);
    out_filler_ = [out](
        not_null<typename Profile::Message*> const message) {
      Profile::Fill(out, message);
    };
  }
}

template<typename Profile>
Method<Profile>::~Method() {
  CHECK(returned_);
  if (Recorder::active_recorder_ != nullptr) {
    serialization::Method method;
    auto* const extension =
        method.MutableExtension(Profile::Message::extension);
    if (out_filler_ != nullptr) {
      out_filler_(extension);
    }
    if (return_filler_ != nullptr) {
      return_filler_(extension);
    }
    Recorder::active_recorder_->WriteAtDestruction(method);
  }
}

template<typename Profile>
template<typename P, typename>
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
    return_filler_ =
        [result](not_null<typename Profile::Message*> const message) {
          Profile::Fill(result, message);
        };
  }
  return result;
}

}  // namespace internal_method
}  // namespace journal
}  // namespace principia
