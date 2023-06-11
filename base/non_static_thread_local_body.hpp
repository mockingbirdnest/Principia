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
    : initialize_([this, args...] {
        return make_not_null_unique<T>(args...);
      }) {}

template<typename T>
template<typename U, typename>
non_static_thread_local<T>::non_static_thread_local(
    std::initializer_list<U> initializer_list)
    : initialize_([this, initializer_list] {
        return make_not_null_unique<T>(initializer_list);
      }) {}

template<typename T>
non_static_thread_local<T>::~non_static_thread_local() {
  absl::MutexLock l1(&RegisteredInstanceMap::global_lock_);
  for (not_null const thread_instances : RegisteredInstanceMap::extant_maps_) {
    absl::MutexLock l2(&thread_instances->thread_lock_);
    auto const mutable_map =
        std::make_shared<InstanceMap>(*thread_instances->map_);
    mutable_map->erase(this);
    thread_instances->map_ = mutable_map;
  }
}

template<typename T>
T const& non_static_thread_local<T>::operator()() const& {
  return Get();
}

template<typename T>
T& non_static_thread_local<T>::operator()() & {
  return Get();
}

template<typename T>
T&& non_static_thread_local<T>::operator()() && {
  return std::move(Get());
}

template<typename T>
non_static_thread_local<T>::RegisteredInstanceMap::RegisteredInstanceMap() {
  absl::MutexLock l(&global_lock_);
  extant_maps_.insert(this);
}

template<typename T>
non_static_thread_local<T>::RegisteredInstanceMap::~RegisteredInstanceMap() {
  absl::MutexLock l(&global_lock_);
  extant_maps_.erase(this);
}

template<typename T>
T& non_static_thread_local<T>::Get() {
  if (!this_thread_all_instances_.has_value()) {
    this_thread_all_instances_.emplace();
  }
  std::shared_ptr<InstanceMap const> const map =
      this_thread_all_instances_->map_;
  // Since we are not locking |map_lock_|, some of the keys in |*map| may point
  // to objects that have since been destroyed, but since we have a consistent
  // snapshot and |*this| cannot have been destroyed, this call to |find| should
  // be fine.
  auto const it = map->find(this);
  if (it != map->end()) {
    return *it->second;
  } else {
    // |thread_lock_| being a member of a thread-local member, this line cannot
    // contend with itself.  It can only contend with the destruction of a
    // |non_static_thread_local<T>|.
    absl::MutexLock l(&this_thread_all_instances_->thread_lock_);
    auto const map =
        std::make_shared<InstanceMap>(*this_thread_all_instances_->map_);
    auto const [it, _] = map->emplace(this, initialize_());
    this_thread_all_instances_->map_ = map;
    return *it->second;
  }
}

template<typename T>
absl::Mutex non_static_thread_local<T>::RegisteredInstanceMap::global_lock_;

template<typename T>
std::set<not_null<typename non_static_thread_local<T>::RegisteredInstanceMap*>>
    non_static_thread_local<T>::RegisteredInstanceMap::extant_maps_;

template<typename T>
thread_local std::optional<
    typename non_static_thread_local<T>::RegisteredInstanceMap>
    non_static_thread_local<T>::this_thread_all_instances_;

}  // namespace internal
}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia
