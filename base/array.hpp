
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
               std::is_convertible<OtherElement*, Element*>::value>::type>
  Array(Array<OtherElement> const& other);
  // No allocation of memory.
  template<typename Size,
           typename =
               typename std::enable_if<std::is_integral<Size>::value>::type>
  Array(Element* data, Size size);

  // Implicit conversion from strings, vectors, and the like.
  template<
      typename Container,
      typename = std::enable_if_t<
          std::is_convertible<decltype(std::declval<Container>().data()),
                              Element*>::value &&
          std::is_integral<decltype(std ::declval<Container>().size())>::value>>
  constexpr Array(Container& container);

  // Construction from a string literal if |Element| is a character type or some
  // flavour of byte.
  template<std::size_t size_plus_1,
           typename Character,
           typename = std::enable_if_t<
               size_plus_1 >= 1 &&
               (std::is_same<Element, unsigned char const>::value ||
                std::is_same<Element, char const>::value ||
                std::is_same<Element, wchar_t const>::value ||
                std::is_same<Element, char16_t const>::value ||
                std::is_same<Element, char32_t const>::value) &&
               (std::is_same<Element, Character>::value ||
                (sizeof(Element) == 1 &&
                 std::is_same<Character, char const>::value))>>
  constexpr explicit Array(Character (&characters)[size_plus_1]);

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

// Deep comparisons.
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral<LeftElement>::value &&
                                     std::is_integral<RightElement>::value>>
bool operator==(Array<LeftElement> left, Array<RightElement> right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral<LeftElement>::value &&
                                     std::is_integral<RightElement>::value>>
bool operator==(Array<LeftElement> left,
                UniqueArray<RightElement> const& right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral<LeftElement>::value &&
                                     std::is_integral<RightElement>::value>>
bool operator==(UniqueArray<LeftElement> const& left,
                Array<RightElement> right);
template<typename LeftElement,
         typename RightElement,
         typename = std::enable_if_t<std::is_integral<LeftElement>::value &&
                                     std::is_integral<RightElement>::value>>
bool operator==(UniqueArray<LeftElement> const& left,
                UniqueArray<RightElement> const& right);

}  // namespace base
}  // namespace principia

#include "base/array_body.hpp"
