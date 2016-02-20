
#include "journal/player.hpp"

#include <string>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "journal/profiles.hpp"
#include "glog/logging.h"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;

namespace journal {

Player::Player(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::in) {
  CHECK(!stream_.fail());
}

bool Player::Play() {
  std::unique_ptr<serialization::Method> method = Read();
  if (method == nullptr) {
    return false;
  }
  // TODO(phl): We don't want to run this method, it directs the output to
  // stderr.log.  Remove it from the protocol buffer at some point.  This
  // will be incompatible with existing journals.
  if (method->HasExtension(serialization::InitGoogleLogging::extension)) {
    return true;
  }

#include "journal/player.generated.cc"
  last_method_ = std::move(method);
  return true;
}

serialization::Method const & Player::last_method() const {
  return *last_method_;
}

std::unique_ptr<serialization::Method> Player::Read() {
  std::string const line = GetLine(&stream_);
  if (line.empty()) {
    return nullptr;
  }

  uint8_t const* const hexadecimal =
      reinterpret_cast<uint8_t const*>(line.c_str());
  int const hexadecimal_size = strlen(line.c_str());
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(),
                               static_cast<int>(bytes.size)));

  return method;
}

}  // namespace journal
}  // namespace principia
