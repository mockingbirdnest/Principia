
#pragma once

#include <cstdint>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

#include "base/array.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "gipfeli/compression.h"
#include "google/protobuf/message.h"
#include "google/protobuf/io/zero_copy_stream.h"

namespace principia {
namespace base {
namespace internal_pull_serializer {

using google::compression::Compressor;

// An output stream based on an array that delegates to a function the handling
// of the case where one array is full.  It calls the |on_full| function passed
// at construction and proceeds with filling the array returned by that
// function.
class DelegatingArrayOutputStream
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  // The stream is supported by |bytes.data| which has size |bytes.size|.  Once
  // that array has been filled, |on_full| is called to somehow consume the
  // data.  |on_full| also returns another array where the stream may output
  // more data.
  DelegatingArrayOutputStream(Bytes bytes,
                              std::function<Bytes(Bytes bytes)> on_full);

  // The ZeroCopyOutputStream API.
  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  Bytes bytes_;
  std::function<Bytes(Bytes bytes)> on_full_;

  std::int64_t byte_count_;
  std::int64_t position_;
  std::int64_t last_returned_size_;  // How many bytes we returned last time
                                     // Next() was called.
};

// This class support serialization which is "pulled" by the client.  That is,
// the client creates a |PullSerializer| object, calls |Start| to start the
// serialization process, repeatedly calls |Pull| to obtain a chunk of data, and
// finally destroys the |PullSerializer|.  |PullSerializer| is intended for use
// in memory-critical contexts as it bounds the amount of memory used
// irrespective of the size of the message to serialize.
class PullSerializer final {
 public:
  // The |size| of the data objects returned by |Pull| are never greater than
  // |chunk_size|.  At most |number_of_chunks| chunks are held in the internal
  // queue.  This class uses at most
  // |number_of_chunks * (chunk_size + O(1)) + O(1)| bytes.
  //TODO(phl):comment
  PullSerializer(int chunk_size, int number_of_chunks, Compressor* compressor);
  ~PullSerializer();

  // Starts the serializer, which will proceed to serialize |message|.  This
  // method must be called at most once for each serializer object.
  void Start(
      not_null<std::unique_ptr<google::protobuf::Message const>> message);

  // Obtain the next chunk of data from the serializer.  Blocks if no data is
  // available.  Returns a |Bytes| object of |size| 0 at the end of the
  // serialization.  The returned object may become invalid the next time |Pull|
  // is called.
  Bytes Pull();

 private:
  // Enqueues the chunk of data to be returned to |Pull| and returns a free
  // chunk.  Blocks if there are no free chunks.  Used as a callback for the
  // underlying |DelegatingArrayOutputStream|.
  Bytes Push(Bytes bytes);

  std::unique_ptr<google::protobuf::Message const> message_;

  Compressor* const compressor_;
  int const chunk_size_;
  int const compressed_chunk_size_;
  int const number_of_chunks_;
  int const number_of_compression_chunks_;

  // The array supporting the stream and the stream itself.
  std::unique_ptr<std::uint8_t[]> data_;
  DelegatingArrayOutputStream stream_;

  // The thread doing the actual serialization.
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |queue_|.
  std::mutex lock_;
  std::condition_variable queue_has_room_;
  std::condition_variable queue_has_elements_;

  // The |queue_| contains the |Bytes| objects filled by |Push| and not yet
  // consumed by |Pull|.  If a |Bytes| object has been handed over to the caller
  // by |Pull| it stays in the queue until the next call to |Pull|, to make sure
  // that the pointer is not reused while the caller processes it.
  std::queue<Bytes> queue_ GUARDED_BY(lock_);

  // The |free_| queue contains the start addresses of chunks that are not yet
  // ready to be returned by |Pull|.  That includes the chunk currently being
  // filled by the stream.
  std::queue<not_null<std::uint8_t*>> free_ GUARDED_BY(lock_);
};

}  // namespace internal_pull_serializer

using internal_pull_serializer::PullSerializer;

}  // namespace base
}  // namespace principia

#include "base/pull_serializer_body.hpp"
