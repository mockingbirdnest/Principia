#pragma once

#include "base/array.hpp"

#include <concepts>
#include <memory>
#include <utility>

#include "glog/logging.h"

namespace principia {
namespace base {
namespace _array {
namespace internal {

template<typename Element>
template<typename OtherElement>
  requires(std::is_convertible_v<OtherElement*, Element*>)
Array<Element>::Array(Array<OtherElement> const& other)
    : data(other.data), size(other.size) {}

template<typename Element>
template<std::integral Size>
Array<Element>::Array(Element* const data, Size const size)
    : data(data), size(static_cast<std::int64_t>(size)) {}

template<typename Element>
template<typename Container>
  requires(std::is_convertible_v<decltype(std::declval<Container>().data()),
                                 Element*> &&
           std::is_integral_v<decltype(std ::declval<Container>().size())>)
constexpr Array<Element>::Array(Container& container)
    : data(container.data()),
      size(static_cast<std::int64_t>(container.size())) {}

template<typename Element>
template<typename Container>
  requires(std::is_convertible_v<
               decltype(std::declval<Container const>().data()),
               Element*> &&
           std::is_integral_v<
               decltype(std ::declval<Container const>().size())>)
constexpr Array<Element>::Array(Container const& container)
    : data(container.data()),
      size(static_cast<std::int64_t>(container.size())) {}

template<typename Element>
template<std::size_t size_plus_1, typename Character>
  requires(size_plus_1 >= 1 &&
           (std::is_same_v<Element, unsigned char const> ||
            std::is_same_v<Element, char const> ||
            std::is_same_v<Element, wchar_t const> ||
            std::is_same_v<Element, char16_t const> ||
            std::is_same_v<Element, char32_t const>) &&
           (std::is_same_v<Element, Character> ||
            (sizeof(Element) == 1 && std::is_same_v<Character, char const>)))
constexpr Array<Element>::Array(Character (&characters)[size_plus_1])
    : data(reinterpret_cast<Element*>(characters)),
      size(size_plus_1 - 1) {
  // The `requires` clause should prevent this from failing, but we explicitly
  // check that the cast is trivial or reinterprets a `char const*`.  The cast
  // is C-style rather than a reinterpret so that this constructor is constexpr
  // in the trivial case.
  static_assert(std::is_same_v<Element, Character> ||
                    std::is_same_v<Character, char const>,
                "reinterpret_cast is unsafe");
  if (characters[size] != 0) {
    LOG(FATAL) << "Array constructed with a character array terminated by the "
                  "non-null 0x"
               << std::hex << static_cast<uint32_t>(characters[size]);
  }
}

template<typename Element>
UniqueArray<Element>::UniqueArray() : size(0) {}

template<typename Element>
template<std::integral Size>
UniqueArray<Element>::UniqueArray(Size const size)
    : data(size == 0 ? nullptr : new Element[static_cast<std::size_t>(size)]),
      size(static_cast<std::int64_t>(size)) {}

template<typename Element>
template<std::integral Size>
UniqueArray<Element>::UniqueArray(std::unique_ptr<Element[]> data,
                                  Size const size)
    : data(data.release()),
      size(static_cast<std::int64_t>(size)) {}

template<typename Element>
Array<Element> UniqueArray<Element>::get() const {
  return Array<Element>(data.get(), size);
}

template<typename Element, std::int32_t max_size>
template<typename... Args>
constexpr BoundedArray<Element, max_size>::BoundedArray(Args&&... args)
    : data_{{std::forward<Args>(args)...}},
      size_(sizeof...(args)) {}

template<typename Element, std::int32_t max_size>
void BoundedArray<Element, max_size>::push_back(const Element& value) {
  data_[size_++] = value;
}

template<typename Element, std::int32_t max_size>
void BoundedArray<Element, max_size>::push_back(Element&& value) {
  data_[size_++] = value;
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::iterator
BoundedArray<Element, max_size>::begin() {
  return data_.begin();
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::iterator
BoundedArray<Element, max_size>::end() {
  return data_.begin() + size_;
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_iterator
BoundedArray<Element, max_size>::begin() const {
  return data_.begin();
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_iterator
BoundedArray<Element, max_size>::end() const {
  return data_.begin() + size_;
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_reverse_iterator
BoundedArray<Element, max_size>::rbegin() const {
  return data_.rend() - size_;
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_reverse_iterator
BoundedArray<Element, max_size>::rend() const {
  return data_.rend();
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::reference
BoundedArray<Element, max_size>::front() {
  return data_.front();
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_reference
BoundedArray<Element, max_size>::front() const {
  return data_.front();
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::reference
BoundedArray<Element, max_size>::back() {
  return data_[size_ - 1];
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::const_reference
BoundedArray<Element, max_size>::back() const {
  return data_[size_ - 1];
}

template<typename Element, std::int32_t max_size>
bool BoundedArray<Element, max_size>::empty() const {
  return size_ == 0;
}

template<typename Element, std::int32_t max_size>
typename BoundedArray<Element, max_size>::size_type
BoundedArray<Element, max_size>::size() const {
  return size_;
}

template<std::integral LeftElement, std::integral RightElement>
bool operator==(Array<LeftElement> left, Array<RightElement> right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<std::size_t>(right.size)) == 0;
}

template<std::integral LeftElement, std::integral RightElement>
bool operator==(Array<LeftElement> left,
                UniqueArray<RightElement> const& right) {
  return left == right.get();
}

template<std::integral LeftElement, std::integral RightElement>
bool operator==(UniqueArray<LeftElement> const& left,
                Array<RightElement> right) {
  return left.get() == right;
}

template<std::integral LeftElement, std::integral RightElement>
bool operator==(UniqueArray<LeftElement> const& left,
                UniqueArray<RightElement> const& right) {
  return left.get() == right.get();
}

}  // namespace internal
}  // namespace _array
}  // namespace base
}  // namespace principia
