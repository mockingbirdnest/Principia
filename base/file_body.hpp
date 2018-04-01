#pragma once

#include "base/file.hpp"

#include <filesystem>
#include <string>

#include "base/macros.hpp"
#include "glog/logging.h"

// In 15.7.0 preview 2 <filesystem> is not yet up to speed compared to
// <experimental/filesystem>.
#if PRINCIPIA_COMPILER_MSVC && _MSC_FULL_VER <= 191426316
#include <experimental/filesystem>
#endif

namespace principia {
namespace base {
namespace internal_file {

inline OFStream::OFStream() {}

inline OFStream::OFStream(std::filesystem::path const& path) {
#if PRINCIPIA_COMPILER_MSVC
  CHECK(path.has_filename()) << path;
  std::filesystem::path directory = path;
  directory.remove_filename();
#if _MSC_FULL_VER <= 191426316
  if (!std::experimental::filesystem::exists(directory.native())) {
    CHECK(std::experimental::filesystem::create_directories(directory.native()))
        << directory;
  }
#else
  if (!std::filesystem::exists(directory)) {
    CHECK(std::filesystem::create_directories(directory))
        << directory;
  }
#endif
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
