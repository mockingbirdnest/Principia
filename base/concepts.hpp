#pragma once

#include <concepts>

#include "serialization/geometry.pb.h"

namespace principia {
namespace base {
namespace _concepts {
namespace internal {

// TODO(phl): This concept applies to a frame.  Extend it.
template<typename F>
concept serializable = requires(serialization::Frame const& message) {
  F::ReadFromMessage(message);
};

}  // namespace internal

using internal::serializable;

}  // namespace _concepts
}  // namespace base
}  // namespace principia
