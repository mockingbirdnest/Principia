#include "journal/recorder.hpp"

#include <filesystem>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "base/serialization.hpp"
#include "base/version.hpp"
#include "glog/logging.h"

namespace principia {
namespace journal {
namespace _recorder {
namespace internal {

using namespace principia::base::_array;
using namespace principia::base::_hexadecimal;
using namespace principia::base::_serialization;
using namespace principia::base::_version;

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

void Recorder::Activate(not_null<Recorder*> const recorder) {
  CHECK(active_recorder_ == nullptr);
  active_recorder_ = recorder;

  // When the recorder gets activated, pretend that we got a GetVersion call.
  // This will record the version at the beginning of the journal, which is
  // useful for forensics.
  serialization::Method method;
  not_null<serialization::GetVersion*> const get_version =
      method.MutableExtension(serialization::GetVersion::extension);
  active_recorder_->WriteAtConstruction(method);
  not_null<serialization::GetVersion::Out*> const out =
      get_version->mutable_out();
  out->set_build_date(BuildDate);
  out->set_version(Version);
  active_recorder_->WriteAtDestruction(method);
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

Recorder* Recorder::active_recorder_ = nullptr;

}  // namespace internal
}  // namespace _recorder
}  // namespace journal
}  // namespace principia
