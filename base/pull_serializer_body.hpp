#pragma once

#include "base/pull_serializer.hpp"

#include <algorithm>

namespace principia {

using std::placeholders::_1;
using std::swap;

namespace base {

namespace internal {

inline DelegatingArrayOutputStream::DelegatingArrayOutputStream(
    Bytes const bytes,
    std::function<Bytes(Bytes const bytes)> on_full)
    : size_(bytes.size),
      data_(bytes.data),
      on_full_(std::move(on_full)),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayOutputStream::Next(void** data, int* size) {
  if (position_ == size_) {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the other array.
    Bytes const bytes = on_full_(Bytes(data_, size_));
    position_ = 0;
    last_returned_size_ = 0;
    data_ = bytes.data;
    size_ = bytes.size;
  }
  CHECK_LT(position_, size_);
  last_returned_size_ = size_ - position_;
  *data = &data_[position_];
  *size = last_returned_size_;
  position_ += last_returned_size_;
  return true;
}

inline void DelegatingArrayOutputStream::BackUp(int count) {
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
    Bytes const bytes = on_full_(Bytes(data_, position_));
    position_ = 0;
    data_ = bytes.data;
    size_ = bytes.size;
  }
  last_returned_size_ = 0;
}

inline std::int64_t DelegatingArrayOutputStream::ByteCount() const {
  return position_;
}

}  // namespace internal

inline PullSerializer::PullSerializer(int const chunk_size,
                                      int const number_of_chunks)
    : chunk_size_(chunk_size),
      number_of_chunks_(number_of_chunks),
      data_(std::make_unique<std::uint8_t[]>(chunk_size_ * number_of_chunks_)),
      stream_(Bytes(data_.get(), chunk_size_),
              std::bind(&PullSerializer::Push, this, _1)) {
  // Mark all the chunks as free except for the first one.
  for (int i = 1; i < number_of_chunks_; ++i) {
    free_.insert(data_.get() + i * chunk_size_);
  }
}

inline PullSerializer::~PullSerializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

inline void PullSerializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {
  CHECK(thread_ == nullptr);
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->SerializeToZeroCopyStream(&stream_));
    // Put a sentinel at the end of the serialized stream so that the client
    // knows that this is the end.  The sentinel doesn't count towards the queue
    // size.
    {
      std::unique_lock<std::mutex> l(lock_);
      queue_.push(Bytes());
    }
    queue_has_elements_.notify_all();
  });
}

inline Bytes PullSerializer::Pull() {
  Bytes result;
  {
    std::unique_lock<std::mutex> l(lock_);
    LOG(ERROR)<<"pulling "<<queue_.size();
    queue_has_elements_.wait(l, [this]() { return !queue_.empty(); });
    //LOG(ERROR)<<"waited "<<queue_.size();
    result = queue_.front();
    queue_.pop();
    free_.insert(result.data);
  }
  LOG(ERROR)<<"pull "<<result.size<<" "<<(void*)(&*result.data);
  queue_has_room_.notify_all();
  return result;
}

inline Bytes PullSerializer::Push(Bytes const bytes) {
  Bytes result;
  CHECK_GE(chunk_size_, bytes.size);
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_room_.wait(l, [this]() {
      return queue_.size() < static_cast<size_t>(number_of_chunks_);
    });
    queue_.emplace(bytes.data, bytes.size);
    CHECK(!free_.empty());
    result = Bytes(free_.front(), chunk_size_);
    free_.pop_front();
    LOG(ERROR)<<"push "<<queue_.size()<<" "<<queue_.back().size;
  }
  queue_has_elements_.notify_all();
  return result;
}

}  // namespace base
}  // namespace principia
