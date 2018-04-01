
#include "journal/recorder.hpp"

#include <filesystem>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "glog/logging.h"

namespace principia {

using base::HexadecimalEncode;
using base::UniqueBytes;

namespace journal {

Recorder::Recorder(std::filesystem::path const& path)
    : stream_(path, std::ios::out) {
  CHECK(!stream_.fail()) << path;
}

Recorder::~Recorder() {
  stream_.close();
}

void Recorder::Write(serialization::Method const& method) {
  CHECK_LT(0, method.ByteSize()) << method.DebugString();
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

void Recorder::Activate(base::not_null<Recorder*> const journal) {
  CHECK(active_recorder_ == nullptr);
  active_recorder_ = journal;
}

void Recorder::Deactivate() {
  CHECK(active_recorder_ != nullptr);
  delete active_recorder_;
  active_recorder_ = nullptr;
}

bool Recorder::IsActivated() {
  return active_recorder_ != nullptr;
}

thread_local Recorder* Recorder::active_recorder_ = nullptr;

}  // namespace journal
}  // namespace principia
