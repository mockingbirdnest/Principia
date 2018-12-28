
#pragma once

#include "base/arena_allocator.hpp"

#include <cstdint>

namespace principia {
namespace base {
namespace arena_allocator_internal {

template<class T>
ArenaAllocator<T>::ArenaAllocator(
    not_null<google::protobuf::Arena*> const arena)
    : arena_(arena) {}

template<class T>
T* ArenaAllocator<T>::allocate(std::size_t const n) {
  return reinterpret_cast<T*>(
      google::protobuf::Arena::CreateArray<std::uint8_t>(arena_,
                                                         n * sizeof(T)));
}

template<class T>
void ArenaAllocator<T>::deallocate(T* const p, std::size_t const n) noexcept {
  // No deallocation, the storage lasts for as long as |arena_|.
}

template<class T, class U>
bool operator==(ArenaAllocator<T> const& at, ArenaAllocator<U> const& au) {
  return at.arena_ == au.arena_;
}

template<class T, class U>
bool operator!=(ArenaAllocator<T> const& at, ArenaAllocator<U> const& au) {
  return at.arena_ != au.arena_;
}

}  // namespace arena_allocator_internal
}  // namespace base
}  // namespace principia
