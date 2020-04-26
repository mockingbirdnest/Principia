#pragma once

#include <type_traits>

namespace principia {
namespace base {
namespace internal_traits {

template<typename... Ts>
constexpr bool all_different_v = false;

template<>
constexpr bool all_different_v<> = true;

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


template<typename, typename = void, typename = void>
struct has_sfinae_read_from_message : std::false_type {};

template<typename T>
struct has_sfinae_read_from_message<
    T, std::void_t<decltype(&T::template ReadFromMessage<>)>>
    : std::true_type {};

template<typename, typename = void, typename = void>
struct has_unconditional_read_from_message : std::false_type {};

template<typename T>
struct has_unconditional_read_from_message<
    T, std::void_t<decltype(&T::ReadFromMessage)>>
    : std::true_type {};

}  // namespace internal_traits

// True if and only if U is an instance of T.
template<template<typename...> typename T, typename U>
inline constexpr bool is_instance_of_v =
    internal_traits::is_instance_of<T, U>::value;

// True if and only if T and U are the same template.
template<template<typename...> typename T, template<typename...> typename U>
inline constexpr bool is_same_template_v =
    internal_traits::is_same_template<T, U>::value;

// True if and only if T has a (possibly templated) static member function named
// ReadFromMessage.
template<typename T>
inline constexpr bool is_serializable_v =
    internal_traits::has_sfinae_read_from_message<T>::value ||
    internal_traits::has_unconditional_read_from_message<T>::value;

using internal_traits::all_different_v;

}  // namespace base
}  // namespace principia
