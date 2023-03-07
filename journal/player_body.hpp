#pragma once

#include "journal/player.hpp"

#include <list>

#include "glog/logging.h"

namespace principia {
namespace journal {
namespace _player {
namespace internal {

template<typename Profile>
bool Player::RunIfAppropriate(serialization::Method const& method_in,
                              serialization::Method const& method_out_return) {
  if (method_in.HasExtension(Profile::Message::extension)) {
    CHECK(method_out_return.HasExtension(Profile::Message::extension))
        << "Unpaired methods:\n"
        << method_in.DebugString() << "\n"
        << method_out_return.DebugString();
    serialization::Method merged_method = method_in;
    merged_method.MergeFrom(method_out_return);
    Profile::Run(merged_method.GetExtension(Profile::Message::extension),
                 pointer_map_);
    return true;
  }
  return false;
}

}  // namespace internal
}  // namespace _player
}  // namespace journal
}  // namespace principia
