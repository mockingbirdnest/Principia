#pragma once

#include "base/file.hpp"

#include <filesystem>
#include <string>
#include <system_error>

#include "base/macros.hpp"
#include "glog/logging.h"

namespace principia {
namespace base {
namespace internal_file {

inline OFStream::OFStream() {}

inline OFStream::OFStream(std::filesystem::path const& path) {
#if PRINCIPIA_COMPILER_MSVC
  CHECK(path.has_filename()) << path;
  std::filesystem::path directory = path;
  directory.remove_filename();
  if (!std::filesystem::exists(directory)) {
    // VS 2017 15.8 Preview 2 has a bug where it returns false if the path ends
    // with a \.
#if _MSC_FULL_VER <= 191526608
    auto d = directory.native();
    d = d.substr(0, d.size() - 1);
    directory = d;
#endif
    std::error_code e;
    CHECK(std::filesystem::create_directories(directory, e))
        << directory << " " << e << " " << e.message();
  }
#endif
  stream_.open(path);
  CHECK(stream_.good()) << path;
}

inline OFStream::~OFStream() {
  stream_.close();
}

inline OFStream& OFStream::operator=(OFStream&& other) {
  CHECK(other.stream_.good());
  stream_ = std::move(other.stream_);
  return *this;
}

inline OFStream&  OFStream::operator<<(std::string const& s) {
  CHECK(stream_.good());
  stream_ << s;
  return *this;
}

}  // namespace internal_file
}  // namespace base
}  // namespace principia
