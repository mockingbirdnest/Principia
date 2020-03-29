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

inline OFStream::OFStream(std::filesystem::path const& path) {
#if PRINCIPIA_COMPILER_MSVC
  CHECK(path.has_filename()) << path;
  // Don't use |remove_filename| here as it leaves a trailing \.  See
  // https://developercommunity.visualstudio.com/content/problem/278829/stdfilesystemcreate-directories-returns-false-if-p.html
  // for a discussion of what happens in that case.
  std::filesystem::path const directory = path.parent_path();
  if (!std::filesystem::exists(directory)) {
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
