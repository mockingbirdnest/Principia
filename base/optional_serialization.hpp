
#pragma once

#include <experimental/optional>

namespace principia {
namespace base {

// We take a MessagePtr rather than a Message* or a not_null<Message*> in order
// to allow either to be deduced.
template<typename T, typename MessagePtr>
void WriteOptionalToMessage(MessagePtr message,
                            std::experimental::optional<T> value);

}  // namespace base
}  // namespace principia

#include "base/optional_serialization_body.hpp"
