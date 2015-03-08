#pragma once

#include "base/push_deserializer.hpp"

#include "glog/logging.h"

namespace principia {

using std::swap;

namespace base {

namespace internal {

inline DelegatingArrayInputStream::DelegatingArrayInputStream(
    std::function<Bytes()> on_empty)
    : on_empty_(std::move(on_empty)),
      byte_count_(0),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayInputStream::Next(void const** const data,
                                             int* const size) {
  //LOG(ERROR)<<"next "<<position_<<" "<<size_;
  if (position_ == bytes_.size) {
    // We're at the end of the array.  Obtain a new one.
    bytes_ = on_empty_();
    //LOG(ERROR)<<"got "<<bytes.size;
    position_ = 0;
    last_returned_size_ = 0;  // Don't let caller back up.
    if (bytes_.size == 0) {
      // At end of input data.
      return false;
    }
  }
  CHECK_LT(position_, bytes_.size);
  last_returned_size_ = bytes_.size - position_;
  *data = &bytes_.data[position_];
  *size = static_cast<int>(last_returned_size_);
  byte_count_ += last_returned_size_;
  position_ += last_returned_size_;
  //LOG(ERROR)<<"ret "<<*size;
  return true;
}

inline void DelegatingArrayInputStream::BackUp(int const count) {
  //LOG(ERROR)<<"backup";
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  byte_count_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up.
}

inline bool DelegatingArrayInputStream::Skip(int const count) {
  //LOG(ERROR)<<"skip";
  CHECK_GE(count, 0);
  last_returned_size_ = 0;   // Don't let caller back up.
  std::int64_t remaining = count;
  while (remaining > bytes_.size - position_) {
    byte_count_ += bytes_.size - position_;
    remaining -= bytes_.size - position_;
    // We're at the end of the array.  Obtain a new one.
    bytes_ = on_empty_();
    position_ = 0;
    if (bytes_.size == 0) {
      // At end of input data.
      return false;
    }
  }
  byte_count_ += remaining;
  position_ += remaining;
  return true;
}

inline std::int64_t DelegatingArrayInputStream::ByteCount() const {
  //LOG(ERROR)<<"bc";
  return byte_count_;
}

}  // namespace internal

inline PushDeserializer::PushDeserializer(int const chunk_size,
                                          int const number_of_chunks)
    : chunk_size_(chunk_size),
      number_of_chunks_(number_of_chunks),
      stream_(std::bind(&PushDeserializer::Pull, this)) {}

inline PushDeserializer::~PushDeserializer() {
  if (thread_ != nullptr) {
    //LOG(ERROR)<<"join";
    thread_->join();
  }
  //LOG(ERROR)<<"destroyed";
}

inline void PushDeserializer::Start(
    not_null<google::protobuf::Message*> const message) {
  CHECK(thread_ == nullptr);
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->ParseFromZeroCopyStream(&stream_));
    //LOG(ERROR)<<"stop";
  });
}

inline void PushDeserializer::Push(Bytes const bytes) {
  // Slice the incoming data in chunks of size at most |chunk_size|.  Release
  // the lock after each chunk to give the deserializer a chance to run.  This
  // method can be called with |bytes| of size 0 to terminate the
  // deserialization, but it never generates a chunk of size 0 in other
  // circumstances.
  Bytes current = bytes;
  CHECK_LE(0, bytes.size);
  do {
    {
      std::unique_lock<std::mutex> l(lock_);
      queue_has_room_.wait(l, [this]() {
        return queue_.size() < static_cast<size_t>(number_of_chunks_);
      });
      queue_.emplace(current.data,
                     std::min(current.size, 
                              static_cast<std::int64_t>(chunk_size_)));
      LOG(ERROR)<<"push "<<queue_.size()<<" "<<queue_.back().size<<" "<<(void*)(&*current.data);
    }
    queue_has_elements_.notify_all();
    current.data = &current.data[chunk_size_];
    current.size -= chunk_size_;
  } while (current.size > 0);
}

inline Bytes PushDeserializer::Pull() {
  Bytes result;
  {
    std::unique_lock<std::mutex> l(lock_);
    LOG(ERROR)<<"pulling "<<queue_.size();
    queue_has_elements_.wait(l, [this]() { return !queue_.empty(); });
    //LOG(ERROR)<<"waited "<<queue_.size();
    result = queue_.front();
    queue_.pop();
  }
  LOG(ERROR)<<"pull "<<result.size<<" "<<(void*)(&*result.data);
  queue_has_room_.notify_all();
  return result;
}

}  // namespace base
}  // namespace principia
