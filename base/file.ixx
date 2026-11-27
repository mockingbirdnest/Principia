export module principia.base.file;

import std;

import "absl/log/check.h";
import "absl/log/log.h";
import "base/macros.hpp";  // 🧙 For PRINCIPIA_COMPILER_MSVC.

import principia.base.file;

namespace principia {
namespace base {
namespace _file {
namespace internal {
export {

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

}  // export
}  // namespace internal

export using internal::OFStream;

}  // namespace _file
}  // namespace base
}  // namespace principia

#include "base/file_body.hpp"
