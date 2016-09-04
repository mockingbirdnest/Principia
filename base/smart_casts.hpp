#pragma once

#include "base/not_null.hpp"

namespace principia {
namespace base {

template<typename P>
struct is_not_null_non_owner : std::false_type {};

template<typename T>
struct is_not_null_non_owner<not_null<T*>> : std::true_type {};

template<typename Result,
         typename Pointer,
         typename = std::enable_if_t<is_not_null_non_owner<Pointer>>>
Result dynamic_cast_not_null(Pointer const pointer);

}  // namespace base
}  // namespace principia
