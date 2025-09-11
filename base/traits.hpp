#pragma once

#include <type_traits>

namespace principia {
namespace base {
namespace _traits {
namespace internal {

template<typename... Ts>
constexpr bool all_different_v = false;

template<>
inline constexpr bool all_different_v<> = true;

template<typename T0, typename... Ts>
constexpr bool all_different_v<T0, Ts...> =
    (!std::is_same_v<T0, Ts> && ...) && all_different_v<Ts...>;


// Has a true `value` member iff the template `T` can be instantiated with the
// arguments `Args...`.
template<template<typename...> typename T, typename... Args>
struct can_be_instantiated {
  template<template<typename...> typename TT, typename = void>
  struct test : std::false_type {};

  template<template<typename...> typename TT>
  struct test<TT, std::void_t<TT<Args...>>> : std::true_type {};

  static constexpr bool value = test<T>::value;
};


// `matches_instantiation<U, can_instantiate, T, Args...>` has a true `value`
// member iff `U` is the same type as `T<Args...>`.  The parameter
// `can_be_instantiated` must indicate whether `T<Args...>` is a valid
// instantiation.
template<typename U,
         bool can_be_instantiated,
         template<typename...> typename T,
         typename... Args>
struct matches_instantiation : std::false_type {};

template<typename U, template<typename...> typename T, typename... Args>
struct matches_instantiation<U, true, T, Args...> {
  static constexpr bool value = std::is_same_v<T<Args...>, U>;
};


// This trait works even if T is an alias template.
template<template<typename...> typename T, typename U>
struct is_instance_of : std::false_type {};

template<template<typename...> typename T,
         template<typename...> typename U,
         typename... Args>
struct is_instance_of<T, U<Args...>>
    : matches_instantiation<U<Args...>,
                            can_be_instantiated<T, Args...>::value,
                            T,
                            Args...> {};


template<template<typename...> typename T,
  template<typename...> typename U>
struct is_same_template : std::false_type {};

template<template<typename...> typename T>
struct is_same_template<T, T> : std::true_type {};


template<typename T, typename T1, typename T2>
struct other_type;

template<typename T>
struct other_type<T, T, T> {
  using type = T;
};

template<typename T1, typename T2>
struct other_type<T1, T1, T2> {
  using type = T2;
};

template<typename T1, typename T2>
struct other_type<T2, T1, T2> {
  using type = T1;
};

template<typename... Ts>
struct tail;

template<typename T>
struct tail<T> {
  using type = T;
  static constexpr T value(T t) {
    return t;
  }
};

template<typename T, typename... Ts>
struct tail<T, Ts...> {
  using type = typename tail<Ts...>::type;
  static constexpr type value(T, Ts... ts) {
    return tail<Ts...>::value(ts...);
  }
};

}  // namespace internal

using internal::all_different_v;
using internal::tail;

// True if and only if U is an instance of T.
template<template<typename...> typename T, typename U>
inline constexpr bool is_instance_of_v = internal::is_instance_of<T, U>::value;

// Note that the order of template parameters is backward from is_instance_of_v.
// This makes it possible to write
//   requires { { expression } -> instance<T>; }
// or
//   template<instance<T> U>
// but is_instance_of_v<T, U> should be preferred in boolean expressions.
// An exception is when defining concepts; there we use instance<U, T> so as to
// get concept-specific error messages.
template<typename U, template<typename...> typename T>
concept instance = is_instance_of_v<T, U>;

// True if and only if T and U are the same template.
template<template<typename...> typename T, template<typename...> typename U>
inline constexpr bool is_same_template_v =
    internal::is_same_template<T, U>::value;

// If T is T1, returns T2.  If T is T2, returns T1.  Otherwise fails.
template<typename T, typename T1, typename T2>
using other_type_t = typename internal::other_type<T, T1, T2>::type;

template<typename... Ts>
using tail_t = typename internal::tail<Ts...>::type;

}  // namespace _traits
}  // namespace base
}  // namespace principia
