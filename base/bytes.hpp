#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.  The constructors don't allocate memory.
template<typename Element>
struct Array {
  Array();  // An object of size 0.
  template<typename OtherElement,
           typename = typename std::enable_if<
               std::is_convertible<OtherElement, Element>::value>::type>
  Array(Array<OtherElement> const& other);
  template<typename Size,
           typename = std::enable_if_t<std::is_integral<Size>::value, Size>>
  Array(Element* const data, Size const size);

  Element* data;
  std::int64_t size;  // In number of elements.
};

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is owned.  The constructors perform allocation.
template<typename Element>
struct UniqueArray {
  UniqueArray();  // An object of size 0.
  template<typename Size,
           typename = std::enable_if_t<std::is_integral<Size>::value, Size>>
  explicit UniqueArray(Size const size);
  template<typename Size,
           typename = std::enable_if_t<std::is_integral<Size>::value, Size>>
  UniqueArray(std::unique_ptr<Element[]> data, Size const size);

  UniqueArray(UniqueArray&& other);
  UniqueArray& operator=(UniqueArray&& other);

  Array<Element> get() const;

  std::unique_ptr<Element[]> data;
  std::int64_t size;  // In number of elements.
};

// Specializations.
using Bytes = Array<std::uint8_t>;
using UniqueBytes = UniqueArray<std::uint8_t>;

// Perform a deep comparison.
template<typename Element>
bool operator==(Array<Element> left, Array<Element> right);
template<typename Element>
bool operator==(Array<Element> left, UniqueArray<Element> const& right);
template<typename Element>
bool operator==(UniqueArray<Element> const& left, Array<Element> right);
template<typename Element>
bool operator==(UniqueArray<Element> const& left,
                UniqueArray<Element> const& right);

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
