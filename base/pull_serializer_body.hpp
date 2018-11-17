
#pragma once

#include "base/pull_serializer.hpp"

#include <algorithm>

#include "base/sink_source.hpp"

namespace principia {
namespace base {
namespace internal_pull_serializer {

using std::placeholders::_1;
using std::swap;

inline DelegatingArrayOutputStream::DelegatingArrayOutputStream(
    Array<std::uint8_t> const bytes,
    std::function<Array<std::uint8_t>(Array<std::uint8_t> const bytes)> on_full)
    : bytes_(bytes),
      on_full_(std::move(on_full)),
      byte_count_(0),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayOutputStream::Next(void** const data,
                                              int* const size) {
  if (position_ == bytes_.size) {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the next one.
    bytes_ = on_full_(bytes_);
    position_ = 0;
  }
  CHECK_LT(position_, bytes_.size);
  last_returned_size_ = bytes_.size - position_;
  *data = &bytes_.data[position_];
  *size = static_cast<int>(last_returned_size_);
  byte_count_ += last_returned_size_;
  position_ += last_returned_size_;
  return true;
}

inline void DelegatingArrayOutputStream::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  byte_count_ -= count;
  position_ -= count;
  // This is called at the end of the stream, in which case we must notify the
  // client about any data remaining in the stream.  If this is called at other
  // times, well, notifying the client doesn't hurt as long as we don't pass a
  // size of 0.
  if (position_ > 0) {
    bytes_ = on_full_(Array<std::uint8_t>(bytes_.data, position_));
    position_ = 0;
  }
  last_returned_size_ = 0;
}

inline std::int64_t DelegatingArrayOutputStream::ByteCount() const {
  return byte_count_;
}

inline PullSerializer::PullSerializer(int const chunk_size,
                                      int const number_of_chunks,
                                      std::unique_ptr<Compressor> compressor)
    : compressor_(std::move(compressor)),
      chunk_size_(chunk_size),
      compressed_chunk_size_(
          compressor_ == nullptr
              ? chunk_size_
              : compressor_->MaxCompressedLength(chunk_size_)),
      number_of_chunks_(number_of_chunks),
      number_of_compression_chunks_(compressor_ == nullptr ? 0 : 1),
      data_(std::make_unique<std::uint8_t[]>(compressed_chunk_size_ *
                                             number_of_chunks_)),
      stream_(Array<std::uint8_t>(data_.get(), chunk_size_),
              std::bind(&PullSerializer::Push, this, _1)) {
  // Check the compatibility of the wait conditions in Push and Pull.
  CHECK_GT(number_of_chunks_ - number_of_compression_chunks_ - 1, 1);

  // Mark all the chunks as free except the last one which is a sentinel for the
  // |queue_|.  The 0th chunk has been passed to the stream, but it's still free
  // until the first call to |on_full|.  Note that the last
  // |compressed_chunk_size_ - chunk_size_| bytes of each chunk are not
  // considered as free.
  for (int i = 0; i < number_of_chunks_ - 1; ++i) {
    free_.push(data_.get() + i * compressed_chunk_size_);
  }
  queue_.push(Array<std::uint8_t>(
      data_.get() + (number_of_chunks_ - 1) * compressed_chunk_size_, 0));
}

inline PullSerializer::~PullSerializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

inline void PullSerializer::Start(
    not_null<std::unique_ptr<google::protobuf::Message const>> message) {
  owned_message_ = std::move(message);
  Start(owned_message_.get());
}

inline void PullSerializer::Start(
    not_null<google::protobuf::Message const*> const message) {
  CHECK(thread_ == nullptr);
  message_ = message;
  thread_ = std::make_unique<std::thread>([this](){
    CHECK(message_->SerializeToZeroCopyStream(&stream_));
    // Put a sentinel at the end of the serialized stream so that the client
    // knows that this is the end.
    Array<std::uint8_t> bytes;
    {
      absl::MutexLock l(&lock_);
      CHECK(!free_.empty());
      bytes = Array<std::uint8_t>(free_.front(), 0);
    }
    Push(bytes);
  });
}

inline Array<std::uint8_t> PullSerializer::Pull() {
  Array<std::uint8_t> result;
  {
    absl::MutexLock l(&lock_);

    // The element at the front of the queue is the one that was last returned
    // by |Pull| and must be dropped and freed.
    auto const queue_has_elements = [this]() { return queue_.size() > 1; };
    lock_.Await(absl::Condition(&queue_has_elements));

    CHECK_LE(2, queue_.size());
    free_.push(queue_.front().data);
    queue_.pop();
    result = queue_.front();
    CHECK_EQ(number_of_chunks_, queue_.size() + free_.size());
  }
  return result;
}

inline Array<std::uint8_t> PullSerializer::Push(Array<std::uint8_t> bytes) {
  Array<std::uint8_t> result;
  CHECK_GE(chunk_size_, bytes.size);
  if (bytes.size > 0 && compressor_ != nullptr) {
    Array<std::uint8_t> compressed_bytes;
    {
      absl::MutexLock l(&lock_);
      CHECK_LE(1 + number_of_compression_chunks_, free_.size());
      free_.pop();
      compressed_bytes =
          Array<std::uint8_t>(free_.front(), compressed_chunk_size_);
      free_.push(bytes.data);
    }
    // We maintain the invariant that the chunk being filled is at the front of
    // the |free_| queue.
    ArraySource<std::uint8_t> source(bytes);
    ArraySink<std::uint8_t> sink(compressed_bytes);
    compressor_->CompressStream(&source, &sink);
    {
      absl::MutexLock l(&lock_);
      bytes = sink.array();
    }
  }
  {
    absl::MutexLock l(&lock_);

    auto const queue_has_room = [this]() {
      // -1 here is because we want to ensure that there is an entry in the
      // free list, in addition to |result| and to
      // |number_of_compression_chunks_| (if present).
      return queue_.size() < static_cast<std::size_t>(number_of_chunks_) -
                                 number_of_compression_chunks_ - 1;
    };
    lock_.Await(absl::Condition(&queue_has_room));

    queue_.emplace(bytes.data, bytes.size);
    CHECK_LE(2 + number_of_compression_chunks_, free_.size());
    CHECK_EQ(free_.front(), bytes.data);
    free_.pop();
    result = Array<std::uint8_t>(free_.front(), chunk_size_);
    CHECK_EQ(number_of_chunks_, queue_.size() + free_.size());
  }
  return result;
}

}  // namespace internal_pull_serializer
}  // namespace base
}  // namespace principia
