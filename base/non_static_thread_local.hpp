#pragma once

#include <functional>
#include <map>
#include <set>
#include <type_traits>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"

namespace principia {
namespace base {
namespace _non_static_thread_local {
namespace internal {

using namespace principia::base::_not_null;

// A non-static data member of type non_static_thread_local<T> functions as a
// non-static thread-local data member of type T: it is an independent object
// for each instance of the enclosing class and for each accessing thread.
// Access is done through the operator().
// For each T, a global lock is taken as part of the following operations:
// 1. For each object of type non_static_thread_local<T>, its destruction;
// 2. For each thread:
//    a. the first access to an object of type non_static_thread_local<T>,
//    b. if a. happened, the destruction of the thread.
// In addition, for each (object, thread) pair, a lock is taken at the first
// access; however, this lock can only be contended by the destruction of an
// object of type non_static_thread_local<T>.
// Access is otherwise lock-free.
template<typename T>
class non_static_thread_local final {
  using InstanceMap = std::map<not_null<non_static_thread_local*>,
                               not_null<std::shared_ptr<T>>>;
 public:
  template<typename... Args>
  non_static_thread_local(Args&&... args);  // NOLINT(runtime/explicit)
  template<typename U,
           typename = std::enable_if_t<
               std::is_constructible_v<T, std::initializer_list<U>>>>
  non_static_thread_local(std::initializer_list<U> initializer_list);
  ~non_static_thread_local();

  T const& operator()() const&;
  T& operator()() &;
  T&& operator()() &&;

 private:
  class RegisteredInstanceMap final {
   public:
    RegisteredInstanceMap();
    ~RegisteredInstanceMap();

    // This map is read lock-free, but mutated under |thread_lock_| for
    // consistency.
    absl::Mutex thread_lock_;
    std::shared_ptr<InstanceMap const> map_ PT_GUARDED_BY(thread_lock_) =
        std::make_shared<InstanceMap const>();

    static absl::Mutex global_lock_ ACQUIRED_BEFORE(thread_lock_);
    static std::set<not_null<RegisteredInstanceMap*>> extant_maps_
        GUARDED_BY(global_lock_);
  };

  T& Get();

  std::function<not_null<std::unique_ptr<T>>()> const initialize_;
  // Emplaced lazily on the first call to |Get| from this thread.
  static thread_local std::optional<RegisteredInstanceMap>
      this_thread_all_instances_;
};

}  // namespace internal

using internal::non_static_thread_local;

}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia

#include "base/non_static_thread_local_body.hpp"
