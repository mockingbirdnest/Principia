
#pragma once

#include "base/optional_serialization.hpp"

namespace principia {
namespace base {

template<typename T, typename Message>
void WriteToOptional(OptionalField<Message> field,
                     std::experimental::optional<T> const& value) {
  if (value) {
    value->WriteToMessage(field.mutable_field());
  } else {
    field.clear_field();
  }
}

}  // namespace base
}  // namespace principia
