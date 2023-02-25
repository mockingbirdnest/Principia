#pragma once

#include <functional>
#include <type_traits>

namespace principia {
namespace base {
namespace _optional_serialization {
namespace internal {

#define SET_OPTIONAL(message_expression, value_identifier)        \
  [&value = (value_identifier), message = (message_expression)] { \
    if (value) {                                                  \
      message->set_##value_identifier(*value);                    \
    } else {                                                      \
      message->clear_##value_identifier();                        \
    }                                                             \
  }()

#define WRITE_TO_OPTIONAL(message_expression, value_identifier)     \
  [&value = (value_identifier), message = (message_expression)]  {  \
    if (value) {                                                    \
      value->WriteToMessage(message->mutable_##value_identifier()); \
    } else {                                                        \
      message->clear_##value_identifier();                          \
    }                                                               \
  }()

}  // namespace internal


}  // namespace _optional_serialization
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_optional_serialization;
}  // namespace principia::base
