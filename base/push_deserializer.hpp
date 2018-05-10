
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
namespace internal_push_deserializer {

using google::compression::Compressor;

// An input stream based on an array that delegates to a function the handling
// of the case where the array is empty.  It calls the |on_empty| function
// passed at construction and proceeds with deserializing the array returned by
// that function.
class DelegatingArrayInputStream
    : public google::protobuf::io::ZeroCopyInputStream {
 public:
  explicit DelegatingArrayInputStream(
      std::function<Array<std::uint8_t>()> on_empty);

  // The ZeroCopyInputStream API.
  bool Next(void const** data, int* size) override;
  void BackUp(int count) override;
  bool Skip(int count) override;
  std::int64_t ByteCount() const override;

 private:
  Array<std::uint8_t> bytes_;
  std::function<Array<std::uint8_t>()> on_empty_;

  std::int64_t byte_count_;
  std::int64_t position_;
  std::int64_t last_returned_size_;  // How many bytes we returned last time
                                     // Next() was called.
};

// This class support deserialization which is "pushed" by the client.  That is,
// the client creates a |PushDeserializer| object, calls |Start| to start the
// deserialization process, repeatedly calls |Push| to send chunks of data for
// deserialization, and finally destroys the |PushDeserializer|.
// |PushDeserializer| is intended for use in memory-critical contexts as it
// bounds the amount of memory used irrespective of the size of the message to
// deserialize.
class PushDeserializer final {
 public:
  // The |size| of the data chunks sent to |Pull| are never greater than
  // |chunk_size|.  The internal queue holds at most |number_of_chunks| chunks.
  // Therefore, this class uses at most
  // |number_of_chunks * (chunk_size + O(1)) + O(1)| bytes.
  PushDeserializer(int chunk_size,
                   int number_of_chunks,
                   std::unique_ptr<Compressor> compressor);
  ~PushDeserializer();

  // Starts the deserializer, which will proceed to deserialize data into
  // |message|.  This method must be called at most once for each deserializer
  // object.  The |done| callback is called once deserialization has completed
  // (which only happens once the client has called |Push| with a chunk of size
  // 0).
  void Start(not_null<std::unique_ptr<google::protobuf::Message>> message,
             std::function<void(google::protobuf::Message const&)> done);

  // Pushes in the internal queue chunks of data that will be extracted by
  // |Pull|.  Splits |bytes| into chunks of at most |chunk_size|.  May block to
  // stay within the maximum size of the queue.  The caller must push an object
  // of size 0 to signal the end of input.  |done| is called once the
  // serialization of |bytes| is complete.  The client must ensure that |bytes|
  // remains live until the call to |done|.  It may reclaim any memory
  // associated with |bytes| in |done|.
  void Push(Array<std::uint8_t> bytes, std::function<void()> done);

  // Same as above but ownership is taken.
  void Push(UniqueArray<std::uint8_t> bytes);

 private:
  // Obtains the next chunk of data from the internal queue.  Blocks if no data
  // is available.  Used as a callback for the underlying
  // |DelegatingArrayOutputStream|.
  Array<std::uint8_t> Pull();

  std::unique_ptr<google::protobuf::Message> message_;

  std::unique_ptr<Compressor> const compressor_;

  // The chunk size passed at construction.  The stream consumes chunks of that
  // size.
  int const chunk_size_;

  // The maximum size of a chunk after compression.  Greater than |chunk_size_|
  // because the compressor will occasionally expand data.  This is the
  // maximum size of the chunks passed to |Push| by the client.
  int const compressed_chunk_size_;

  // The number of chunks passed at construction, used to size |data_|.
  int const number_of_chunks_;

  UniqueArray<std::uint8_t> uncompressed_data_;

  DelegatingArrayInputStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |queue_|.
  std::mutex lock_;
  std::condition_variable queue_has_room_;
  std::condition_variable queue_has_elements_;
  // The |queue_| contains the |Array<std::uint8_t>| object filled by |Push| and
  // not yet consumed by |Pull|.  The |done_| queue contains the callbacks.  The
  // two queues are out of step: an element is removed from |queue_| by |Pull|
  // when it returns a chunk to the stream, but the corresponding callback is
  // removed from |done_| (and executed) when |Pull| returns.
  std::queue<Array<std::uint8_t>> queue_ GUARDED_BY(lock_);
  std::queue<std::function<void()>> done_ GUARDED_BY(lock_);
};

}  // namespace internal_push_deserializer

using internal_push_deserializer::PushDeserializer;

}  // namespace base
}  // namespace principia

#include "base/push_deserializer_body.hpp"
