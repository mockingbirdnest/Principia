
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

// A simple container for an array and a size.
template<typename Element, std::int32_t max_size>
class BoundedArray final {
  using Container = std::array<Element, max_size>;

 public:
  using iterator = typename Container::iterator;
  using const_iterator = typename Container::const_iterator;
  using const_reverse_iterator = typename  Container::const_reverse_iterator;
  using reference = typename Container::reference;
  using const_reference = typename Container::const_reference;
  using size_type = typename Container::size_type;
  using value_type = Element;

  template<typename... Args>
  constexpr BoundedArray(Args&&... args);

  void push_back(const Element& value);
  void push_back(Element&& value);

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

  reference front();
  const_reference front() const;
  reference back();
  const_reference back() const;

  bool empty() const;
  size_type size() const;

 private:
  Container data_;
  std::int32_t size_;
};


// Specializations.
using Bytes = Array<std::uint8_t>;
using UniqueBytes = UniqueArray<std::uint8_t>;

// Swap.
template<typename Element>
void swap(Array<Element>& left, Array<Element>& right);

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
