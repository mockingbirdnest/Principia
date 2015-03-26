#pragma once

#include "base/bytes.hpp"

#include "glog/logging.h"

namespace principia {
namespace base {

template<typename T>
UniqueArray<T>::UniqueArray() : size(0) {}

template<typename T>
template<typename U, typename>
UniqueArray<T>::UniqueArray(U const size)
    : data(new std::uint8_t[static_cast<size_t>(size)]),
      size(static_cast<std::int64_t>(size)) {}

template<typename T>
template<typename U, typename>
UniqueArray<T>::UniqueArray(std::unique_ptr<T[]> data, U const size)
    : data(data.release()),
      size(static_cast<std::int64_t>(size)) {}

template<typename T>
UniqueArray<T>::UniqueArray(UniqueArray&& other)
    : data(std::move(other.data)), size(other.size) {}

template<typename T>
UniqueArray<T>& UniqueArray<T>::operator=(UniqueArray&& other) {
  data = std::move(other.data);
  size = other.size;
  return *this;
}

template<typename T>
Array<T> UniqueArray<T>::get() const {
  return Array<T>(data.get(), size);
}

template<typename T>
Array<T>::Array() : data(nullptr), size(0) {}

template<typename T>
template<typename W, typename>
Array<T>::Array(Array<W> const& other) : data(other.data), size(other.size) {}

template<typename T>
template<typename U, typename>
Array<T>::Array(T* const data, U const size)
    : data(data), size(static_cast<std::int64_t>(size)) {}

template<typename T>
bool operator==(Array<T> left, Array<T> right) {
  if (left.size != right.size) {
    return false;
  }
  return std::memcmp(left.data,
                     right.data,
                     static_cast<size_t>(right.size)) == 0;
}

template<typename T>
inline bool operator==(Array<T> left, UniqueArray<T> const& right) {
  return left == right.get();
}

template<typename T>
inline bool operator==(UniqueArray<T> const& left, Array<T> right) {
  return left.get() == right;
}

template<typename T>
inline bool operator==(UniqueArray<T> const& left, UniqueArray<T> const& right) {
  return left.get() == right.get();
}

}  // namespace base
}  // namespace principia
