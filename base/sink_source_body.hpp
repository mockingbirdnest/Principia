#pragma once

#include "base/sink_source.hpp"

#include <algorithm>

namespace principia {
namespace base {
namespace _sink_source {
namespace internal {

template<typename Element>
ArraySource<Element>::ArraySource(Array<Element> const& array)
    : array_(array) {}

template<typename Element>
std::size_t ArraySource<Element>::Available() const {
  return array_.size - next_to_read_;
}

template<typename Element>
const char* ArraySource<Element>::Peek(std::size_t* const length) {
  *length = array_.size - next_to_read_;
  return reinterpret_cast<const char*>(array_.data + next_to_read_);
}

template<typename Element>
void ArraySource<Element>::Skip(std::size_t const n) {
  next_to_read_ += n;
}

template<typename Element>
ArraySink<Element>::ArraySink(Array<Element> const& array) : array_(array) {}

template<typename Element>
Array<Element> ArraySink<Element>::array() const {
  Array<Element> result;
  result.data = array_.data;
  result.size = next_to_write_;
  return result;
}

template<typename Element>
void ArraySink<Element>::Append(const char* const data, std::size_t const n) {
  // Do no copying if the caller filled in the result of GetAppendBuffer()
  if (data != reinterpret_cast<const char*>(array_.data + next_to_write_)) {
    memcpy(array_.data + next_to_write_, data, n);
  }
  next_to_write_ += n;
}

template<typename Element>
char* ArraySink<Element>::GetAppendBuffer(
    std::size_t const min_size,
    std::size_t const desired_size_hint,
    char* const scratch,
    std::size_t const scratch_size,
    std::size_t* const allocated_size) {
  *allocated_size = std::min(static_cast<std::int64_t>(desired_size_hint),
                             array_.size - next_to_write_);
  return reinterpret_cast<char*>(array_.data + next_to_write_);
}

}  // namespace internal
}  // namespace _sink_source
}  // namespace base
}  // namespace principia
