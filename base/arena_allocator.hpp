
#pragma once

#include <cstddef>

#include "base/not_null.hpp"
#include "google/protobuf/arena.h"

namespace principia {
namespace base {
namespace arena_allocator_internal {

template<class T>
class ArenaAllocator {
 public:
  using value_type = T;

  explicit ArenaAllocator(not_null<google::protobuf::Arena*> arena);

  template<class U>
  ArenaAllocator(ArenaAllocator<U> const& au) noexcept;

  T* allocate(std::size_t n);

  void deallocate(T* p, std::size_t) noexcept;

 private:
  not_null<google::protobuf::Arena*> const arena_;

  template<class T, class U>
  friend bool operator==(ArenaAllocator<T> const& at,
                         ArenaAllocator<U> const& au);
  template<class T, class U>
  friend bool operator!=(ArenaAllocator<T> const& at,
                         ArenaAllocator<U> const& au);

  template<class U>
  friend class ArenaAllocator;
};

template<class T, class U>
bool operator==(ArenaAllocator<T> const& at, ArenaAllocator<U> const& au);
template<class T, class U>
bool operator!=(ArenaAllocator<T> const& at, ArenaAllocator<U> const& au);

}  // namespace arena_allocator_internal

using arena_allocator_internal::ArenaAllocator;

}  // namespace base
}  // namespace principia

#include "base/arena_allocator_body.hpp"
