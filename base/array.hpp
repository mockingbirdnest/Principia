#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {
namespace _array {
namespace internal {

// A simple container for a pointer and size.  `data` is not owned.
template<typename Element>
struct Array final {
  // An object of size 0.
  Array() = default;
  // Mostly useful for adding constness.
  template<typename OtherElement,
           typename = typename std::enable_if_t<
               std::is_convertible_v<OtherElement*, Element*>>>
  Array(Array<OtherElement> const& other);
  // No allocation of memory.
  template<typename Size,
           typename =
               typename std::enable_if_t<std::is_integral_v<Size>>>
  Array(Element* data, Size size);

  // Implicit conversion from strings, vectors, and the like.
  template<
      typename Container,
      typename = std::enable_if_t<
          std::is_convertible_v<decltype(std::declval<Container>().data()),
                              Element*> &&
          std::is_integral_v<decltype(std ::declval<Container>().size())>>>
  constexpr Array(Container& container);  // NOLINT(runtime/explicit)

  template<typename Container,
           typename = std::enable_if_t<
               std::is_convertible_v<
                   decltype(std::declval<Container const>().data()),
                   Element*> &&
               std::is_integral_v<
                   decltype(std ::declval<Container const>().size())>>>
  constexpr Array(Container const& container);  // NOLINT(runtime/explicit)

  // Construction from a string literal if `Element` is a character type or some
  // flavour of byte.
  template<std::size_t size_plus_1,
           typename Character,
           typename = std::enable_if_t<
               size_plus_1 >= 1 &&
               (std::is_same_v<Element, unsigned char const> ||
                std::is_same_v<Element, char const> ||
                std::is_same_v<Element, wchar_t const> ||
                std::is_same_v<Element, char16_t const> ||
                std::is_same_v<Element, char32_t const>) &&
               (std::is_same_v<Element, Character> ||
                (sizeof(Element) == 1 &&
                 std::is_same_v<Character, char const>))>>
  constexpr explicit Array(Character (&characters)[size_plus_1]);

  Element* data;
  std::int64_t size = 0;  // In number of elements.
};

// A simple container for a pointer and size.  `data` is owned.
template<typename Element>
struct UniqueArray final {
  // An object of size 0.
  UniqueArray();
  // Allocates memory for `size` elements.
  template<typename Size,
           typename =
               typename std::enable_if_t<std::is_integral_v<Size>>>
  explicit UniqueArray(Size size);
  // Takes ownership of an existing array.
  template<typename Size,
           typename =
               typename std::enable_if_t<std::is_integral_v<Size>>>
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
  constexpr BoundedArray(Args&&... args);  // NOLINT(runtime/explicit)

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

// Deep comparisons.
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral_v<LeftElement> &&
                                     std::is_integral_v<RightElement>>>
bool operator==(Array<LeftElement> left, Array<RightElement> right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral_v<LeftElement> &&
                                     std::is_integral_v<RightElement>>>
bool operator==(Array<LeftElement> left,
                UniqueArray<RightElement> const& right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral_v<LeftElement> &&
                                     std::is_integral_v<RightElement>>>
bool operator==(UniqueArray<LeftElement> const& left,
                Array<RightElement> right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral_v<LeftElement> &&
                                     std::is_integral_v<RightElement>>>
bool operator==(UniqueArray<LeftElement> const& left,
                UniqueArray<RightElement> const& right);

}  // namespace internal

using internal::Array;
using internal::BoundedArray;
using internal::operator==;
using internal::UniqueArray;

}  // namespace _array
}  // namespace base
}  // namespace principia

#include "base/array_body.hpp"
