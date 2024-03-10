#pragma once

#include <concepts>

namespace principia {
namespace base {
namespace _concepts {
namespace internal {

template<typename T>
concept serializable = requires {
  &T::ReadFromMessage;
} || requires {
  &T::template ReadFromMessage<>;
};

}  // namespace internal

using internal::serializable;

}  // namespace _concepts
}  // namespace base
}  // namespace principia
