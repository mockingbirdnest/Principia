#pragma once

#include <experimental/filesystem>
#include <fstream>

namespace principia {
namespace base {
namespace internal_file {

// A RAII wrapper for a |std::ofstream|.
class OFStream {
 public:
  OFStream();
  explicit OFStream(std::experimental::filesystem::path const& path);
  ~OFStream();

  OFStream& operator=(OFStream&& other);
  OFStream& operator<<(std::string const& s);

 private:
  std::ofstream stream_;
};

}  // namespace internal_file

using internal_file::OFStream;

}  // namespace base
}  // namespace principia

#include "base/file_body.hpp"
