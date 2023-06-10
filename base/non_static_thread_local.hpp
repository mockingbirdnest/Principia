#pragma once

#include <functional>
#include <map>
#include <set>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"

namespace principia {
namespace base {
namespace _non_static_thread_local {
namespace internal {

using namespace principia::base::_not_null;

template<typename T>
class non_static_thread_local final {
 public:
  template<typename... Args>
  non_static_thread_local(Args&&... args);
  ~non_static_thread_local();

  T const& operator()() const&;
  T& operator()() &;
  T&& operator()() &&;

 private:
  class MemberMap {
   public:
    MemberMap() {
      absl::MutexLock l(&lock_);
      extant_maps_.insert(this);
    }

    ~MemberMap() {
      absl::MutexLock l(&lock_);
      extant_maps_.erase(this);
    }

    std::map<not_null<non_static_thread_local*>, T> map_;

    static absl::Mutex lock_;
    static std::set<not_null<MemberMap*>> extant_maps_ GUARDED_BY(&lock_);
  };
  std::function<T&()> const get_;
  static thread_local MemberMap members_;
};

template<typename T>
absl::Mutex non_static_thread_local<T>::MemberMap::lock_;
template<typename T>
std::set<not_null<typename non_static_thread_local<T>::MemberMap*>>
    non_static_thread_local<T>::MemberMap::extant_maps_;

template<typename T>
thread_local typename non_static_thread_local<T>::MemberMap
    non_static_thread_local<T>::members_;

}  // namespace internal

using internal::non_static_thread_local;

}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia

#include "base/non_static_thread_local_body.hpp"
