#pragma once

#include "base/push_deserializer.hpp"

#include "glog/logging.h"

namespace principia {

using std::swap;

namespace base {

namespace internal {

MyStream::MyStream(not_null<std::uint8_t*> data,
                   int const size,
                   std::function<Bytes()> on_empty)
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
    Bytes const bytes = on_empty_();
    //TODO(phl): and then?
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

PushDeserializer::PushDeserializer(int const max_size)
    : data_(std::make_unique<std::uint8_t[]>(max_size << 1)),
      stream_(data_.get(),
              max_size,
              std::bind(&PushDeserializer::Pull, this)),
      done_(false) {}

PushDeserializer::~PushDeserializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

void PushDeserializer::Start(
    not_null<google::protobuf::Message*> const message) {
  CHECK(thread_ == nullptr);
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->ParseFromZeroCopyStream(&stream_));
    {
      std::unique_lock<std::mutex> l(lock_);
      done_ = true;
    }
    holder_is_full_.notify_all();
  });
}

void PushDeserializer::Push(Bytes const bytes) {
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_empty_.wait(l, [this]() { return holder_ == nullptr; });
    holder_ = std::make_unique<Bytes>(bytes.data, bytes.size);
  }
  holder_is_full_.notify_all();
}

Bytes PushDeserializer::Pull() {
  std::unique_ptr<Bytes const> result;
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_full_.wait(l, [this]() { return done_ || holder_ != nullptr; });
    if (holder_ != nullptr) {
      result = std::move(holder_);
    }
  }
  holder_is_empty_.notify_all();
  if (result == nullptr) {
    // Done.
    return Bytes::Null;
  } else {
    return *result;
  }
}

}  // namespace base
}  // namespace principia
