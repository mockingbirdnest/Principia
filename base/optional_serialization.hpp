
#pragma once

#include <experimental/optional>
#include <functional>
#include <type_traits>

namespace principia {
namespace base {

template<typename Message>
struct OptionalField {
  std::function<Message*()> mutable_field;
  std::function<void()> clear_field;
};

#define OPTIONAL_FIELD(message, field)                      \
  (::principia::base::OptionalField<                        \
       std::remove_cv_t<std::remove_reference_t<            \
           decltype(message->field())>>> {                  \
       [m = (message)]() { return m->mutable_##field(); },  \
       [m = (message)]() { m->clear_##field(); }})

template<typename T, typename Message>
void WriteToOptional(OptionalField<Message> field,
                     std::experimental::optional<T> const& value);

}  // namespace base
}  // namespace principia

#include "base/optional_serialization_body.hpp"
