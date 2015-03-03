#pragma once

#include "base/push_deserializer.hpp"

namespace principia {

using std::swap;

namespace base {

namespace internal {

MyStream::MyStream(not_null<std::uint8_t*> data,
                   int const size,
                   std::function<void(Bytes const bytes)> on_empty)
    : size_(size),
      data1_(&data[0]),
      data2_(&data[size_]),
      on_empty_(std::move(on_empty)),
      position_(0),
      last_returned_size_(0) {}

bool MyStream::Next(const void** data, int* size) {
  if (position_ == size_) {
    // We're at the end of the array.  Hand the current array over to the
    // callback to be filled.
    //TODO(phl): What if it doesn't fill it completely?
    on_empty_(Bytes(data1_, size_));
    position_ = 0;
    last_returned_size_ = 0;
    swap(data1_, data2_);
  }
  CHECK_LT(position_, size_);
  last_returned_size_ = size_ - position_;
  *data = &data1_[position_];
  *size = last_returned_size_;
  position_ += last_returned_size_;
  return true;
}

void MyStream::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up further.
}

bool MyStream::Skip(int count) {
  CHECK_GE(count, 0);
  last_returned_size_ = 0;   // Don't let caller back up.
  if (count > size_ - position_) {
    position_ = size_;
    return false;
  } else {
    position_ += count;
    return true;
  }
}

std::int64_t MyStream::ByteCount() const {
  return position_;
}

}  // namespace internal

}  // namespace base
}  // namespace principia
