#pragma once

#include <iterator>

namespace principia {
namespace base {
namespace _ranges {
namespace internal {

// A range [begin, end[, suitable for use in a range-based for loop.
template<typename Iterator>
class IterableRange {
  using IteratorTraits = std::iterator_traits<Iterator>;

 public:
  static_assert(
      std::is_base_of<std::forward_iterator_tag,
                      typename IteratorTraits::iterator_category>::value,
      "Iterator must be a forward iterator");
  using iterator = Iterator;
  // There is no way to obtain the const_iterator from the iterator; we could
  // have something specialized for all iterators of standard containers, but
  // for now we simply restrict ourselves to immutable ranges.
  static_assert(
      std::is_const<
          std::remove_reference_t<decltype(*std::declval<Iterator>())>>::value,
      "Iterator must be a const iterator");
  using const_iterator = Iterator;
  using difference_type = typename IteratorTraits::difference_type;
  using size_type = std::make_unsigned_t<difference_type>;
  using value_type = typename IteratorTraits::value_type;
  using reference = typename IteratorTraits::reference;
  using const_reference = value_type const&;

  IterableRange(Iterator begin, Iterator end);

  const_iterator begin() const;
  const_iterator end() const;

  const_iterator cbegin() const;
  const_iterator cend() const;

  size_type size() const;
  bool empty() const;

 private:
  Iterator begin_;
  Iterator end_;
};

template<typename Iterator>
IterableRange<Iterator> Range(Iterator begin, Iterator end);

}  // namespace internal

using internal::IterableRange;
using internal::Range;

}  // namespace _ranges
}  // namespace base
}  // namespace principia

#include "base/ranges_body.hpp"
