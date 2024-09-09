#pragma once

#include <filesystem>
#include <fstream>
#include <string>

namespace principia {
namespace base {
namespace _file {
namespace internal {

// A RAII wrapper for a `std::ofstream`.
class OFStream {
 public:
  OFStream() = default;
  explicit OFStream(std::filesystem::path const& path);
  ~OFStream();

  void Flush();

  OFStream& operator=(OFStream&& other);
  OFStream& operator<<(std::string const& s);

 private:
  std::ofstream stream_;
};

}  // namespace internal

using internal::OFStream;

}  // namespace _file
}  // namespace base
}  // namespace principia

#include "base/file_body.hpp"
