#pragma once

#include "base/file.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {
namespace internal_file {

inline OFStream::OFStream() {}

inline OFStream::OFStream(std::experimental::filesystem::path const& path) {
  CHECK(path.has_filename()) << path;
  std::experimental::filesystem::path directory = path;
  directory.remove_filename();
  if (!std::experimental::filesystem::exists(directory)) {
    CHECK(std::experimental::filesystem::create_directories(directory))
        << directory;
  }
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
