#pragma once

#include "base/push_deserializer.hpp"

#include "glog/logging.h"

namespace principia {

using std::swap;

namespace base {

namespace internal {

DelegatingArrayInputStream::DelegatingArrayInputStream(
    std::function<Bytes()> on_empty)
    : size_(0),
      data_(Bytes().data),
      on_empty_(std::move(on_empty)),
      byte_count_(0),
      position_(0),
      last_returned_size_(0) {}

bool DelegatingArrayInputStream::Next(void const** const data,
                                      int* const size) {
  if (position_ == size_) {
    // We're at the end of the array.  Obtain a new one.
    Bytes const bytes = on_empty_();
    size_ = bytes.size;
    data_ = bytes.data;
    position_ = 0;
    last_returned_size_ = 0;  // Don't let caller back up.
    if (size_ == 0) {
      // At end of input data.
      return false;
    }
  }
  CHECK_LT(position_, size_);
  last_returned_size_ = size_ - position_;
  *data = &data_[position_];
  *size = last_returned_size_;
  byte_count_ += last_returned_size_;
  position_ += last_returned_size_;
  return true;
}

void DelegatingArrayInputStream::BackUp(int const count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  byte_count_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up.
}

bool DelegatingArrayInputStream::Skip(int const count) {
  CHECK_GE(count, 0);
  last_returned_size_ = 0;   // Don't let caller back up.
  int remaining = count;
  while (remaining > size_ - position_) {
    byte_count_ += size_ - position_;
    remaining -= size_ - position_;
    // We're at the end of the array.  Obtain a new one.
    Bytes const bytes = on_empty_();
    size_ = bytes.size;
    data_ = bytes.data;
    position_ = 0;
    if (size_ == 0) {
      // At end of input data.
      return false;
    }
  }
  byte_count_ += remaining;
  position_ += remaining;
  return true;
}

std::int64_t DelegatingArrayInputStream::ByteCount() const {
  return byte_count_;
}

}  // namespace internal

PushDeserializer::PushDeserializer(int const chunk_size,
                                   int const number_of_chunks)
    : chunk_size_(chunk_size),
      number_of_chunks_(number_of_chunks),
      stream_(std::bind(&PushDeserializer::Pull, this)) {}

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
      // Append a sentinel.  It doesn't count against the queue capacity.
      std::unique_lock<std::mutex> l(lock_);
      queue_.emplace(Bytes().data, 0);
      //TODO(phl):Does that work?
    }
    queue_has_elements_.notify_all();
  });
}

void PushDeserializer::Push(Bytes const bytes) {
  // Slice the incoming data in chunks of size at most |chunk_size|.  Release
  // the lock after each chunk to give the deserializer a chance to run.
  Bytes current = bytes;
  while (current.size > 0) {
    {
      std::unique_lock<std::mutex> l(lock_);
      queue_has_room_.wait(l, [this]() {
        return queue_.size() < static_cast<size_t>(number_of_chunks_);
      });
      queue_.emplace(current.data, std::min(current.size, chunk_size_));
    }
    queue_has_elements_.notify_all();
    current.data = &current.data[chunk_size_];
    current.size -= chunk_size_;
  }
}

Bytes PushDeserializer::Pull() {
  Bytes result;
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_elements_.wait(l, [this]() { return !queue_.empty(); });
    result = queue_.front();
    queue_.pop();
  }
  queue_has_room_.notify_all();
  return result;
}

}  // namespace base
}  // namespace principia
