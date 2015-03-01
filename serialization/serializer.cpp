#include "serialization/serializer.hpp"

namespace principia {
namespace serialization {

SynchronizingArrayOutputString::SynchronizingArrayOutputString(
    void* data, int size, int block_size)
    : data_(reinterpret_cast<std::uint8_t*>(data)),
      size_(size),
      block_size_(block_size > 0 ? block_size : size),
      position_(0),
      last_returned_size_(0) {}

bool SynchronizingArrayOutputString::Next(void** data, int* size) {
  if (position_ < size_) {
    last_returned_size_ = std::min(block_size_, size_ - position_);
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
  GOOGLE_CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  GOOGLE_CHECK_LE(count, last_returned_size_);
  GOOGLE_CHECK_GE(count, 0);
  position_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up further.
}

std::int64_t SynchronizingArrayOutputString::ByteCount() const {
  return position_;
}

Serializer::Serializer(int const chunk_size) {}

void Serializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {}

base::not_null<std::string const*> Serializer::Get() {}

}  // namespace serialization
}  // namespace principia
