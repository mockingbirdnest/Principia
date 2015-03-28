#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

// A simple container for a pointer and size.  |data| is not owned.
template<typename Element>
struct Array {
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
  Array(Element* const data, Size const size);

  Element* data;
  std::int64_t size;  // In number of elements.
};

// A simple container for a pointer and size.  |data| is owned.
template<typename Element>
struct UniqueArray {
  // An object of size 0.
  UniqueArray();
  // Allocates memory for |size| elements.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  explicit UniqueArray(Size const size);
  // Takes ownership of an existing array.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  UniqueArray(std::unique_ptr<Element[]> data, Size const size);

  // Move it, move it!
  // TODO(phl): The following should be = default when MSVC supports it.
  UniqueArray(UniqueArray&& other);  // NOLINT(build/c++11)
  UniqueArray& operator=(UniqueArray&& other);  // NOLINT(build/c++11)

  // No transfer of ownership.
  Array<Element> get() const;

  std::unique_ptr<Element[]> data;
  std::int64_t size;  // In number of elements.
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
