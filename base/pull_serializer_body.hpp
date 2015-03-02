#pragma once

#include "base/pull_serializer.hpp"

#include <algorithm>

namespace principia {

using std::placeholders::_1;
using std::swap;

namespace base {

namespace internal {

DelegatingTwoArrayOutputStream::DelegatingTwoArrayOutputStream(
    base::not_null<std::uint8_t*> data,
    int const size,
    std::function<void(Bytes const bytes)> on_full)
    : size_(size),
      data1_(&data[0]),
      data2_(&data[size_]),
      on_full_(std::move(on_full)),
      position_(0),
      last_returned_size_(0) {}

bool DelegatingTwoArrayOutputStream::Next(void** data, int* size) {
  if (position_ == size_) {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the other array.
    on_full_(Bytes(data1_, size_));
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

void DelegatingTwoArrayOutputStream::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  // This is called at the end of the stream, in which case we must notify the
  // client about any data remaining in the stream.  If this is called at other
  // times, well, notifying the client doesn't hurt as long as we don't pass a
  // size of 0.
  if (position_ > 0) {
    on_full_(Bytes(data1_, position_));
    position_ = 0;
  }
  last_returned_size_ = 0;
  swap(data1_, data2_);
}

std::int64_t DelegatingTwoArrayOutputStream::ByteCount() const {
  return position_;
}

}  // namespace internal

PullSerializer::PullSerializer(int const max_size)
    : data_(std::make_unique<std::uint8_t[]>(max_size << 1)),
      stream_(data_.get(),
              max_size,
              std::bind(&PullSerializer::Push, this, _1)),
      done_(false) {}

PullSerializer::~PullSerializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

void PullSerializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {
  CHECK(thread_ == nullptr);
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->SerializeToZeroCopyStream(&stream_));
    std::unique_lock<std::mutex> l(lock_);
    done_ = true;
  });
}

Bytes PullSerializer::Pull() {
  std::unique_ptr<Bytes const> result;
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_full_.wait(l, [this](){ return holder_ != nullptr; });
    result = std::move(holder_);
  }
  holder_is_empty_.notify_all();
  return *result;
}

bool PullSerializer::done() const {
  std::unique_lock<std::mutex> l(lock_);
  return done_ && holder_ == nullptr;
}

void PullSerializer::Push(Bytes const bytes) {
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_empty_.wait(l, [this](){ return holder_ == nullptr; });
    holder_ = std::make_unique<Bytes>(bytes.data, bytes.size);
  }
  holder_is_full_.notify_all();
}

}  // namespace base
}  // namespace principia
