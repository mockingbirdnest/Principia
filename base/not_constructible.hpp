#pragma once

namespace principia {
namespace base {
namespace _not_constructible {
namespace internal {

// Type traits should inherit from this, so that they may use inheritance
// without virtual destructors.
struct not_constructible {
  not_constructible() = delete;
  ~not_constructible() = delete;
};

}  // namespace internal

using internal::not_constructible;

}  // namespace _not_constructible
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_not_constructible;
}  // namespace principia::base
