#pragma once

#include <utility>

#include "base/non_static_thread_local.hpp"

namespace principia {
namespace base {
namespace _non_static_thread_local {
namespace internal {

template<typename T>
template<typename... Args>
non_static_thread_local<T>::non_static_thread_local(Args&&... args)
    : get_([this, tuple = std::tuple{std::forward<Args>(args)...}]() -> T& {
        if (!members_.has_value()) {
          members_.emplace();
        }
        auto const [it, _] = members_->map_.emplace(
            std::piecewise_construct, std::tuple{this}, tuple);
        return it->second;
      }) {}

template<typename T>
template<typename U, typename>
non_static_thread_local<T>::non_static_thread_local(
    std::initializer_list<U> initializer_list)
    : get_([this, initializer_list]() -> T& {
        if (!members_.has_value()) {
          members_.emplace();
        }
        auto const [it, _] =
            members_->map_.emplace(std::piecewise_construct,
                                   std::tuple{this},
                                   std::tuple{initializer_list});
        return it->second;
      }) {}

template<typename T>
non_static_thread_local<T>::~non_static_thread_local() {
  absl::MutexLock l(&MemberMap::lock_);
  for (not_null const member_map : MemberMap::extant_maps_) {
    member_map->map_.erase(this);
  }
}

template<typename T>
T const& non_static_thread_local<T>::operator()() const& {
  return get_();
}

template<typename T>
T& non_static_thread_local<T>::operator()() & {
  return get_();
}

template<typename T>
T&& non_static_thread_local<T>::operator()() && {
  return std::move(get_());
}

template<typename T>
non_static_thread_local<T>::MemberMap::MemberMap() {
  absl::MutexLock l(&lock_);
  extant_maps_.insert(this);
}

template<typename T>
non_static_thread_local<T>::MemberMap::~MemberMap() {
  absl::MutexLock l(&lock_);
  extant_maps_.erase(this);
}

template<typename T>
absl::Mutex non_static_thread_local<T>::MemberMap::lock_;

template<typename T>
std::set<not_null<typename non_static_thread_local<T>::MemberMap*>>
    non_static_thread_local<T>::MemberMap::extant_maps_;

template<typename T>
thread_local std::optional<typename non_static_thread_local<T>::MemberMap>
    non_static_thread_local<T>::members_;

}  // namespace internal
}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia
