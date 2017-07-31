
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

// A simple container for a pointer and size.  |data| is not owned.
template<typename Element>
struct Array final {
  // An object of size 0.
  Array();
  // Mostly useful for adding constness.
  template<typename OtherElement,
           typename = typename std::enable_if<
               std::is_convertible<OtherElement, Element>::value>::type>
  Array(Array<OtherElement> const& other);
  // No allocation of memory.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  Array(Element* data, Size size);

  Element* data;
  std::int64_t size;  // In number of elements.
};

// A simple container for a pointer and size.  |data| is owned.
template<typename Element>
struct UniqueArray final {
  // An object of size 0.
  UniqueArray();
  // Allocates memory for |size| elements.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  explicit UniqueArray(Size size);
  // Takes ownership of an existing array.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  UniqueArray(std::unique_ptr<Element[]> data, Size size);

  // Move it, move it!
  UniqueArray(UniqueArray&& other) = default;
  UniqueArray& operator=(UniqueArray&& other) = default;

  // No transfer of ownership.
  Array<Element> get() const;

  std::unique_ptr<Element[]> data;
  std::int64_t size;  // In number of elements.
};

// A simple container for an array and a size.  The client is expected to use
// aggregate initialization for this type, and to ensure that the values passed
// for |data| and |size| are consistent.
template<typename Element, std::int32_t size_>
struct BoundedArray final {
  typename std::array<Element, size_>::const_iterator begin() const;
  typename std::array<Element, size_>::const_iterator end() const;

  std::array<Element, size_> data;
  std::int32_t size;
};


// Specializations.
using Bytes = Array<std::uint8_t>;
using UniqueBytes = UniqueArray<std::uint8_t>;

// Deep comparisons.
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

#include "base/array_body.hpp"
