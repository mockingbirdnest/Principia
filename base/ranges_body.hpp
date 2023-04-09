#pragma once

#include "base/ranges.hpp"

namespace principia {
namespace base {
namespace _ranges {
namespace internal {

template<typename Iterator>
IterableRange<Iterator>::IterableRange(Iterator const begin, Iterator const end)
    : begin_(begin), end_(end) {}

template<typename Iterator>
typename IterableRange<Iterator>::const_iterator
IterableRange<Iterator>::begin() const {
  return begin_;
}

template<typename Iterator>
typename IterableRange<Iterator>::const_iterator IterableRange<Iterator>::end()
    const {
  return end_;
}

template<typename Iterator>
typename IterableRange<Iterator>::const_iterator
IterableRange<Iterator>::cbegin() const {
  return begin_;
}

template<typename Iterator>
typename IterableRange<Iterator>::const_iterator IterableRange<Iterator>::cend()
    const {
  return end_;
}

template<typename Iterator>
typename IterableRange<Iterator>::size_type IterableRange<Iterator>::size()
    const {
  return std::distance(begin_, end_);
}

template<typename Iterator>
bool IterableRange<Iterator>::empty() const {
  return begin_ == end_;
}

template<typename Iterator>
IterableRange<Iterator> Range(Iterator const begin, Iterator const end) {
  return IterableRange<Iterator>(begin, end);
}

}  // namespace internal
}  // namespace _ranges
}  // namespace base
}  // namespace principia
