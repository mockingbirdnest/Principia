#pragma once

namespace principia {
namespace base {
namespace _concepts {
namespace internal {

// True if and only if T has a (possibly templated) static member function named
// ReadFromMessage.
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
