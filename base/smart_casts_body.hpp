#pragma once

#include "base/smart_casts.hpp"

namespace principia {
namespace base {

template<typename Result, typename Pointer, typename>
Result dynamic_cast_not_null(Pointer const pointer) {
  return dynamic_cast<Result>(static_cast<typename Pointer::pointer>(pointer));
}

}  // namespace base
}  // namespace principia
