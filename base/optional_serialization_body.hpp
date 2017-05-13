
#pragma once

#include "base/optional_serialization.hpp"

namespace principia {
namespace base {

template<typename T, typename MessagePtr>
void WriteOptionalToMessage(MessagePtr message,
                            std::experimental::optional<T> value) {
  if (value) {
    value->WriteToMessage(message);
  }
}

}  // namespace base
}  // namespace principia
