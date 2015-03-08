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
    : bytes_(bytes),
      on_full_(std::move(on_full)),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayOutputStream::Next(void** data, int* size) {
  if (position_ == bytes_.size) {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the next one.
    bytes_ = on_full_(bytes_);
    position_ = 0;
    last_returned_size_ = 0;
  }
  CHECK_LT(position_, bytes_.size);
  last_returned_size_ = bytes_.size - position_;
  *data = &bytes_.data[position_];
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
    bytes_ = on_full_(Bytes(bytes_.data, position_));
    position_ = 0;
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
              std::bind(&PullSerializer::Push, this, _1)),
      is_first_pull_(true) {
  // Mark all the chunks as free.  The 0th chunk has been passed to the strem,
  // but it's still free until the first call to |on_full|.
  for (int i = 0; i < number_of_chunks_; ++i) {
    free_.push(data_.get() + i * chunk_size_);
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
    queue_has_elements_.wait(l, [this]() { return !queue_.empty(); });
    // The element at the front of the queue is the one that was last returned
    // by |Pull| and must be dropped and freed, except the first time |Pull| is
    // called.
    if (is_first_pull_) {
      is_first_pull_ = false;
    } else {
      free_.push(queue_.front().data);
      queue_.pop();
    }
    result = queue_.front();
  }
  queue_has_room_.notify_all();
  return result;
}

inline Bytes PullSerializer::Push(Bytes const bytes) {
  Bytes result;
  CHECK_GE(chunk_size_, bytes.size);
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_room_.wait(l, [this]() {
      // -1 here is because we want to ensure that there is an entry in the
      // free list.
      return queue_.size() < static_cast<size_t>(number_of_chunks_) - 1;
    });
    queue_.emplace(bytes.data, bytes.size);
    CHECK(!free_.empty());
    CHECK_EQ(free_.front(), bytes.data);
    free_.pop();
    result = Bytes(free_.front(), chunk_size_);
  }
  queue_has_elements_.notify_all();
  return result;
}

}  // namespace base
}  // namespace principia
