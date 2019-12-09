#pragma once

#include <type_traits>

namespace principia {
namespace base {

// If T has a member named WriteToMessage, provides the member constant |value|
// equal to true. Otherwise |value| is false.
template<typename, typename = std::void_t<>>
struct is_serializable : std::false_type {};

template<typename T>
struct is_serializable<
    T,
    std::void_t<decltype(std::declval<T>().WriteToMessage(nullptr))>>
    : std::true_type {};

template<typename T>
inline constexpr bool is_serializable_v = is_serializable<T>::value;

}  // namespace base
}  // namespace principia
