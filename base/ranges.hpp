#pragma once

namespace principia {
namespace base {

// A range [begin, end[, suitable for use in a range-based for loop.
template<typename Iterator>
class IterableRange {
 public:
  IterableRange(Iterator begin, Iterator end) : begin_(begin), end_(end) {}
  Iterator begin() {
    return begin_;
  }
  Iterator end() {
    return end_;
  }

 private:
  Iterator begin_;
  Iterator end_;
};

template<typename Iterator>
IterableRange<Iterator> Range(Iterator begin, Iterator end) {
  return IterableRange<Iterator>(begin, end);
}

}  // namespace base
}  // namespace principia
