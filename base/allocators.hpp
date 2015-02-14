#pragma once

#include <memory>
#include <type_traits>

namespace principia {
namespace base {

// The following class is adapted from http://stackoverflow.com/a/21028912.

// Allocator adaptor that interposes |construct| calls to convert value
// initialization into default initialization.
// Using this allocator in a container of ints yields a ~10% performance gain in
// calls to |resize()| or construction with a given size.
// http://en.cppreference.com/w/cpp/language/value_initialization,
// http://en.cppreference.com/w/cpp/language/default_initialization,
// http://en.cppreference.com/w/cpp/container/vector/resize,
// http://en.cppreference.com/w/cpp/container/vector/vector,
// http://en.cppreference.com/w/cpp/concept/DefaultInsertable,
// http://en.cppreference.com/w/cpp/concept/MoveInsertable.
template<typename T, typename A = std::allocator<T>>
class DefaultInitializationAllocator : public A {
  using a_t = std::allocator_traits<A>;
 public:

  // NOTE(egg): When MSVC supports it, this class should inherit the
  // constructors of |A|.  In the MSVC implementation of |std::allocator| they
  // do nothing.
  // using A::A;

  template<typename U>
  struct rebind {
    using other = DefaultInitializationAllocator<
                      U, typename a_t::template rebind_alloc<U>>;
  };

  template<typename U>
  void construct(U* ptr) {
    ::new(static_cast<void*>(ptr)) U;
  }

  template<typename U, typename... Args>
  void construct(U* ptr, Args&&... args) {
    a_t::construct(static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
  }
};

}  // namespace base
}  // namespace principia
