#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

template<typename T>
struct Array;//TODO(phl):rorer

//UODO(phl):comment.
template<typename T>
struct UniqueArray {
  UniqueArray();  // An object of size 0.
  template<typename U,
           typename = std::enable_if_t<std::is_integral<U>::value, U>>
  explicit UniqueArray(U const size);
  template<typename U,
           typename = std::enable_if_t<std::is_integral<U>::value, U>>
  UniqueArray(std::unique_ptr<T[]> data, U const size);

  UniqueArray(UniqueArray&& other);
  UniqueArray& operator=(UniqueArray&& other);

  Array<T> get() const;

  std::unique_ptr<T[]> data;
  std::int64_t size;  // In number of |T| elements.
};

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
template<typename T>//UODO(phl): names
struct Array {
  Array();  // An object of size 0.
  template<typename W,
           typename = typename std::enable_if<
               std::is_convertible<W, T>::value>::type>
  Array(Array<W> const& other);
  template<typename U,
           typename = std::enable_if_t<std::is_integral<U>::value, U>>
  Array(T* const data, U const size);

  T* data;
  std::int64_t size;  // In number of |T| elements.
};

// Specializations.
using Bytes = Array<std::uint8_t>;
using UniqueBytes = UniqueArray<std::uint8_t>;

// Performs a deep comparison.
template<typename T>
bool operator==(Array<T> left, Array<T> right);
template<typename T>
bool operator==(Array<T> left, UniqueArray<T> const& right);
template<typename T>
bool operator==(UniqueArray<T> const& left, Array<T> right);
template<typename T>
bool operator==(UniqueArray<T> const& left, UniqueArray<T> const& right);

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
