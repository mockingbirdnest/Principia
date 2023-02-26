#pragma once

namespace principia {
namespace base {
namespace _tags {
namespace internal {

struct uninitialized_t {};
constexpr uninitialized_t uninitialized;

}  // namespace internal

using internal::uninitialized;
using internal::uninitialized_t;

}  // namespace _tags
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_tags;
}  // namespace principia::base
