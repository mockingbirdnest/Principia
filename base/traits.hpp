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


template<template<typename...> typename T, typename U>
struct is_instance_of : std::false_type {};

template<template<typename...> typename T, typename... Args>
struct is_instance_of<T, T<Args...>> : std::true_type {};


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

}  // namespace internal

using internal::all_different_v;

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

}  // namespace _traits
}  // namespace base
}  // namespace principia
