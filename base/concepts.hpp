#pragma once

namespace principia {
namespace base {
namespace _concepts {

// True if and only if T has a (possibly templated) static member function named
// ReadFromMessage.
template<typename T>
concept serializable = requires {
  &T::ReadFromMessage;
} || requires {
  &T::template ReadFromMessage<>;
};

}  // namespace _concepts
}  // namespace base
}  // namespace principia
