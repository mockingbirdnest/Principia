
#pragma once

#include "base/push_deserializer.hpp"

#include <algorithm>

#include "base/sink_source.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/coded_stream_inl.h"

namespace principia {
namespace base {
namespace internal_push_deserializer {

inline DelegatingArrayInputStream::DelegatingArrayInputStream(
    std::function<Bytes()> on_empty)
    : on_empty_(std::move(on_empty)),
      byte_count_(0),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayInputStream::Next(void const** const data,
                                             int* const size) {
  if (position_ == bytes_.size) {
    // We're at the end of the array.  Obtain a new one.
    bytes_ = on_empty_();
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
  return true;
}

inline void DelegatingArrayInputStream::BackUp(int const count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  byte_count_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up.
}

inline bool DelegatingArrayInputStream::Skip(int const count) {
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
  return byte_count_;
}

inline PushDeserializer::PushDeserializer(
    int const chunk_size,
    int const number_of_chunks,
    std::unique_ptr<Compressor> compressor)
    : compressor_(std::move(compressor)),
      chunk_size_(chunk_size),
      compressed_chunk_size_(
          compressor_ == nullptr
              ? chunk_size_
              : compressor_->MaxCompressedLength(chunk_size_)),
      number_of_chunks_(number_of_chunks),
      uncompressed_data_(chunk_size_),
      stream_(std::bind(&PushDeserializer::Pull, this)) {
  // This sentinel ensures that the two queue are correctly out of step.
  done_.push(nullptr);
}

inline PushDeserializer::~PushDeserializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

inline void PushDeserializer::Start(
    not_null<std::unique_ptr<google::protobuf::Message>> message,
    std::function<void(google::protobuf::Message const&)> done) {
  CHECK(thread_ == nullptr);
  message_ = std::move(message);
  thread_ = std::make_unique<std::thread>([this, done]() {
    // It is a well-known annoyance that, in order to set the total byte limit,
    // we have to copy code from MessageLite::ParseFromZeroCopyStream.  Blame
    // Kenton.
    google::protobuf::io::CodedInputStream decoder(&stream_);
    decoder.SetTotalBytesLimit(1 << 29, 1<< 29);
    CHECK(message_->ParseFromCodedStream(&decoder));
    CHECK(decoder.ConsumedEntireMessage());

    // Run any remainining chunk callback.
    std::unique_lock<std::mutex> l(lock_);
    CHECK_EQ(1, done_.size());
    auto const done_front = done_.front();
    if (done_front != nullptr) {
      done_front();
    }
    done_.pop();

    // Run the final callback.
    if (done != nullptr) {
      done(*message_);
    }
  });
}

inline void PushDeserializer::Push(Bytes const bytes,
                                   std::function<void()> done) {
  // Slice the incoming data in chunks of size at most |chunk_size|.  Release
  // the lock after each chunk to give the deserializer a chance to run.  This
  // method should be called with |bytes| of size 0 to terminate the
  // deserialization, but it never generates a chunk of size 0 in other
  // circumstances.  The |done| callback is attached to the last chunk.
  Bytes current = bytes;
  CHECK_LE(0, bytes.size);

  // Decide how much data we are going to push on the queue.  In the presence of
  // compression we have to respect the boundary of the incoming block.  In the
  // absence of compression we have a stream so we can cut into as many chunks
  // as we like.
  int queued_chunk_size;
  if (compressor_ == nullptr) {
    queued_chunk_size = chunk_size_;
  } else {
    CHECK_LE(bytes.size, compressed_chunk_size_);
    queued_chunk_size = compressed_chunk_size_;
  }

  bool is_last;
  do {
    {
      is_last = current.size <= queued_chunk_size;
      std::unique_lock<std::mutex> l(lock_);
      queue_has_room_.wait(l, [this]() {
        return queue_.size() < static_cast<std::size_t>(number_of_chunks_);
      });
      queue_.emplace(current.data,
                     std::min(current.size,
                              static_cast<std::int64_t>(queued_chunk_size)));
      done_.emplace(is_last ? std::move(done) : nullptr);
    }
    queue_has_elements_.notify_all();
    current.data = &current.data[queued_chunk_size];
    current.size -= queued_chunk_size;
  } while (!is_last);
}

inline void PushDeserializer::Push(UniqueBytes bytes) {
  Bytes const unowned_bytes = bytes.get();
  bytes.data.release();
  Push(unowned_bytes,
       /*done=*/[b = unowned_bytes.data]() { delete[] b; });
}

inline Bytes PushDeserializer::Pull() {
  Bytes result;
  {
    std::unique_lock<std::mutex> l(lock_);
    queue_has_elements_.wait(l, [this]() { return !queue_.empty(); });
    // The front of |done_| is the callback for the |Bytes| object that was just
    // processed.  Run it now.
    CHECK(!done_.empty());
    auto const done = done_.front();
    if (done != nullptr) {
      done();
    }
    done_.pop();
    // Get the next |Bytes| object to process and remove it from |queue_|.
    // Uncompress it if needed.
    auto const& front = queue_.front();
    if (front.size == 0 || compressor_ == nullptr) {
      result = front;
    } else {
      ArraySource<std::uint8_t> source(front);
      ArraySink<std::uint8_t> sink(uncompressed_data_.get());
      CHECK(compressor_->UncompressStream(&source, &sink));
      result = sink.array();
    }
    queue_.pop();
  }
  queue_has_room_.notify_all();
  return result;
}

}  // namespace internal_push_deserializer
}  // namespace base
}  // namespace principia
