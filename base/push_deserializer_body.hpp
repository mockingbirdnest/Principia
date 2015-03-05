#pragma once

#include "base/push_deserializer.hpp"

#include "glog/logging.h"

namespace principia {

using std::swap;

namespace base {

namespace internal {

MyStream::MyStream(std::function<Bytes()> on_empty)
    : size_(0),
      data_(Bytes().data),
      on_empty_(std::move(on_empty)),
      position_(0),
      last_returned_size_(0) {}

bool MyStream::Next(const void** data, int* size) {
  if (position_ == size_) {
    // We're at the end of the array.  Obtain a new one.
    Bytes const bytes = on_empty_();
    if (bytes.size == 0) {
      // At end of input data.
      return false;
    }
    size_ = bytes.size;
    data_ = bytes.data;
    position_ = 0;
    last_returned_size_ = 0;
  }
  CHECK_LT(position_, size_);
  last_returned_size_ = size_ - position_;
  *data = &data_[position_];
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

PushDeserializer::PushDeserializer(int const number_of_slots)
    : number_of_slots_(number_of_slots),
      stream_(std::bind(&PushDeserializer::Pull, this)),
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
    queue_has_elements_.notify_all();
  });
}

void PushDeserializer::Push(Bytes const bytes) {
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_room_.wait(l, [this]() {
      return queue_.size() < number_of_slots_;
    });
    queue_.emplace(bytes.data, bytes.size);
  }
  queue_has_elements_.notify_all();
}

Bytes PushDeserializer::Pull() {
  Bytes result;
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_elements_.wait(l, [this]() { return done_ || !queue_.empty(); });

    if (!queue_.empty()) {
      result = queue_.front();
    }
  }
  queue_has_room_.notify_all();
  return result;
}

}  // namespace base
}  // namespace principia
