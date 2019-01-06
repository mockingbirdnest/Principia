
#include "journal/recorder.hpp"

#include <filesystem>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "glog/logging.h"
#include "base/serialization.hpp"

namespace principia {

using base::HexadecimalEncoder;
using base::SerializeAsBytes;
using base::UniqueArray;

namespace journal {

Recorder::Recorder(std::filesystem::path const& path)
    : stream_(path, std::ios::out) {
  CHECK(!stream_.fail()) << path;
}

void Recorder::WriteAtConstruction(serialization::Method const& method) {
  lock_.Lock();
  WriteLocked(method);
}

void Recorder::WriteAtDestruction(serialization::Method const& method) {
  WriteLocked(method);
  lock_.Unlock();
}

void Recorder::Activate(base::not_null<Recorder*> const recorder) {
  CHECK(active_recorder_ == nullptr);
  active_recorder_ = recorder;
}

void Recorder::Deactivate() {
  CHECK(active_recorder_ != nullptr);
  delete active_recorder_;
  active_recorder_ = nullptr;
}

bool Recorder::IsActivated() {
  return active_recorder_ != nullptr;
}

void Recorder::WriteLocked(serialization::Method const& method) {
  static auto* const encoder = new HexadecimalEncoder</*null_terminated=*/true>;
  CHECK_LT(0, method.ByteSize()) << method.DebugString();
  auto const hexadecimal = encoder->Encode(SerializeAsBytes(method).get());
  stream_ << hexadecimal.data.get() << "\n";
  stream_.flush();
}

thread_local Recorder* Recorder::active_recorder_ = nullptr;

}  // namespace journal
}  // namespace principia
