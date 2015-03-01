#pragma once

#include "base/pull_serializer.hpp"

namespace principia {

using std::placeholders::_1;
using std::placeholders::_2;
using std::swap;

namespace base {

DelegatingTwoArrayOutputStream::DelegatingTwoArrayOutputStream(
    base::not_null<std::uint8_t*> data,
    int const size,
    std::function<void(base::not_null<std::uint8_t const*> const data,
                       int const size)> on_full)
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
    on_full_(data1_, size_);
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
  last_returned_size_ = 0;  // Don't let caller back up further.
  //TODO(phl):odd 
  on_full_(data1_, position_);
}

std::int64_t DelegatingTwoArrayOutputStream::ByteCount() const {
  return position_;
}

PullSerializer::Data::Data(base::not_null<std::uint8_t const*> const data,
                           int const size)
    : data(data), size(size) {}

PullSerializer::PullSerializer(int const max_size)
    : data_(std::make_unique<std::uint8_t[]>(max_size << 1)),
      stream_(data_.get(),
              max_size,
              std::bind(&PullSerializer::Push, this, _1, _2)) {}

PullSerializer::~PullSerializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

void PullSerializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {
  CHECK(thread_ == nullptr);
  done_ = false;
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->SerializeToZeroCopyStream(&stream_));
    std::unique_lock<std::mutex> l(lock_);
    done_ = true;
    holder_is_full_.notify_all();
  });
}

PullSerializer::Data PullSerializer::Pull() {
  Data* result;
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_full_.wait(l, [this](){ return done_ || holder_ != nullptr; });
    if (holder_ == nullptr) {
      // Done.
      result = no_data_;
    } else {
      result = holder_.release();
    }
  }
  holder_is_empty_.notify_all();
  return *result;
}

void PullSerializer::Push(base::not_null<std::uint8_t const*> const data,
                          int const size) {
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_empty_.wait(l, [this](){ return holder_ == nullptr; });
    holder_ = std::make_unique<Data>(data, size);
  }
  holder_is_full_.notify_all();
}

std::uint8_t PullSerializer::no_data_data_ = 0;
PullSerializer::Data* PullSerializer::no_data_ =
    new Data(&PullSerializer::no_data_data_, 0);

}  // namespace base
}  // namespace principia
