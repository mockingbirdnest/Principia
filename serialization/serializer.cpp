#include "serialization/serializer.hpp"

namespace principia {
namespace serialization {

SynchronizingArrayOutputString::SynchronizingArrayOutputString(
    std::uint8_t* data, int size)
    : data_(data),
      size_(size),
      position_(0),
      last_returned_size_(0) {}

bool SynchronizingArrayOutputString::Next(void** data, int* size) {
  if (position_ < size_) {
    last_returned_size_ = size_ - position_;
    *data = data_ + position_;
    *size = last_returned_size_;
    position_ += last_returned_size_;
    return true;
  } else {
    // We're at the end of the array.
    last_returned_size_ = 0;   // Don't let caller back up.
    return false;
  }
}

void SynchronizingArrayOutputString::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up further.
}

std::int64_t SynchronizingArrayOutputString::ByteCount() const {
  return position_;
}

Serializer::Serializer(int const chunk_size)
    : data_(std::make_unique<std::uint8_t[]>(chunk_size)),
      stream_(data_.get(), chunk_size) {}

void Serializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {}

base::not_null<std::string const*> Serializer::Get() {}

}  // namespace serialization
}  // namespace principia
