#pragma once

namespace principia {
namespace base {
namespace _tags {
namespace internal {

struct uninitialized_t {};
constexpr uninitialized_t uninitialized;

template<typename T>
concept uninitialized_constructible = requires {
  { T(uninitialized) } -> std::same_as<T>;
};

}  // namespace internal

using internal::uninitialized;
using internal::uninitialized_constructible;
using internal::uninitialized_t;

}  // namespace _tags
}  // namespace base
}  // namespace principia
