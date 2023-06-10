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
// 1. For each object of type non_static_thread_local<T>: destruction;
// 2. For each thread, the first access to an object of type
//    non_static_thread_local<T>.
// Access is otherwise lock-free.
template<typename T>
class non_static_thread_local final {
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
  class MemberMap final {
   public:
    MemberMap();
    ~MemberMap();

    std::map<not_null<non_static_thread_local*>, T> map_;

    static absl::Mutex lock_;
    static std::set<not_null<MemberMap*>> extant_maps_ GUARDED_BY(&lock_);
  };
  std::function<T&()> const get_;
  // Emplaced lazily on the first call to |get_| from this thread.
  static thread_local std::optional<MemberMap> members_;
};

}  // namespace internal

using internal::non_static_thread_local;

}  // namespace _non_static_thread_local
}  // namespace base
}  // namespace principia

#include "base/non_static_thread_local_body.hpp"
