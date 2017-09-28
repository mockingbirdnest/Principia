
#pragma once

#include "base/array.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename Element>
Array<Element>::Array() : data(nullptr), size(0) {}

template<typename Element>
template<typename OtherElement, typename>
Array<Element>::Array(Array<OtherElement> const& other)
    : data(other.data), size(other.size) {}

template<typename Element>
template<typename Size, typename>
Array<Element>::Array(Element* const data, Size const size)
    : data(data), size(static_cast<std::int64_t>(size)) {}

template<typename Element>
UniqueArray<Element>::UniqueArray() : size(0) {}

template<typename Element>
template<typename Size, typename>
UniqueArray<Element>::UniqueArray(Size const size)
    : data(new std::uint8_t[static_cast<std::size_t>(size)]),
      size(static_cast<std::int64_t>(size)) {}

template<typename Element>
template<typename Size, typename>
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

template<typename Element>
bool operator==(Array<Element> left, Array<Element> right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<std::size_t>(right.size)) == 0;
}

template<typename Element>
bool operator==(Array<Element> left, UniqueArray<Element> const& right) {
  return left == right.get();
}

template<typename Element>
bool operator==(UniqueArray<Element> const& left, Array<Element> right) {
  return left.get() == right;
}

template<typename Element>
bool operator==(UniqueArray<Element> const& left,
                UniqueArray<Element> const& right) {
  return left.get() == right.get();
}

}  // namespace base
}  // namespace principia
