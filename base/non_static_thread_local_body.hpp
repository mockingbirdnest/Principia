#pragma once

#include <utility>

#include "base/non_static_thread_local.hpp"

namespace principia {
namespace base {
namespace _non_static_thread_local {
namespace internal {

template<typename T, typename That>
template<typename... Args>
non_static_thread_local<T, That>::non_static_thread_local(
    not_null<That*> const that, Args&&... args)
    : that_(that),
      get_(
          [that, tuple = std::tuple{std::forward<Args>(args)...}]() -> T& {
        auto const [it, inserted] = members_.map_.emplace(
            std::piecewise_construct, std::tuple{that}, tuple);
        return it->second;
      }) {}

template<typename T, typename That>
non_static_thread_local<T, That>::~non_static_thread_local() {
  absl::MutexLock l(&members_.lock_);
  for (not_null const member_map : members_.extant_maps_) {
    member_map->map_.erase(that_);
  }
}

template<typename T, typename That>
T const& non_static_thread_local<T, That>::operator()() const& {
  return get_();
}

template<typename T, typename That>
T& non_static_thread_local<T, That>::operator()() & {
  return get_();
}

template<typename T, typename That>

T&& non_static_thread_local<T, That>::operator()() && {
  return std::move(get_());
}

}  // namespace internal
}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia
