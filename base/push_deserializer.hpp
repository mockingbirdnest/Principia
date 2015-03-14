#pragma once

#include <cstdint>
#include <condition_variable>  // NOLINT(build/c++11)
#include <memory>
#include <mutex>  // NOLINT(build/c++11)
#include <queue>
#include <thread>  // NOLINT(build/c++11)

#include "base/bytes.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "google/protobuf/message.h"
#include "google/protobuf/io/zero_copy_stream.h"

namespace principia {
namespace base {

namespace internal {
// An input stream based on an array that delegates to a function the handling
// of the case where the array is empty.  It calls the |on_empty| function
// passed at construction and proceeds with deserializing the array returned by
// that function.
class DelegatingArrayInputStream
    : public google::protobuf::io::ZeroCopyInputStream {
 public:
  explicit DelegatingArrayInputStream(std::function<Bytes()> on_empty);
  ~DelegatingArrayInputStream() = default;

  // The ZeroCopyInputStream API.
  bool Next(void const** data, int* size) override;
  void BackUp(int count) override;
  bool Skip(int count) override;
  std::int64_t ByteCount() const override;

 private:
  Bytes bytes_;
  std::function<Bytes()> on_empty_;

  std::int64_t byte_count_;
  std::int64_t position_;
  std::int64_t last_returned_size_;  // How many bytes we returned last time
                                     // Next() was called.
};

}  // namespace internal

// This class support deserialization which is "pushed" by the client.  That is,
// the client creates a |PushDeserializer| object, calls |Start| to start the
// deserialization process, repeatedly calls |Push| to send chunks of data for
// deserialization, and finally destroys the |PushDeserializer|.
// |PushDeserializer| is intended for use in memory-critical contexts as it
// bounds the amount of memory used irrespective of the size of the message to
// deserialize.
class PushDeserializer {
 public:
  // The |size| of the data chunks sent to |Pull| are never greater than
  // |chunk_size|.  The internal queue holds at most |number_of_chunks| chunks.
  // Therefore, this class uses at most
  // |number_of_chunks * (chunk_size + O(1)) + O(1)| bytes.
  PushDeserializer(int const chunk_size,
                   int const number_of_chunks);
  ~PushDeserializer();

  // Starts the deserializer, which will proceed to deserialize data into
  // |message|.  This method must be called at most once for each deserializer
  // object.
  void Start(not_null<google::protobuf::Message*> const message);

  // Pushes in the internal queue chunks of data that will be extracted by
  // |Pull|.  Splits |bytes| into chunks of at most |chunk_size|.  May block to
  // stay within the maximum size of the queue.  The caller must push an object
  // of size 0 to signal the end of input.  |done| is called once the
  // serialization of |bytes| is complete.  The client must ensure that |bytes|
  // remains live until the call to |done|.  It may reclaim any memory
  // associated with |bytes| in |done|.
  void Push(Bytes const bytes, std::function<void()> done);

 private:
  // Obtains the next chunk of data from the internal queue.  Blocks if no data
  // is available.  Used as a callback for the underlying
  // |DelegatingArrayOutputStream|.
  Bytes Pull();

  int const chunk_size_;
  int const number_of_chunks_;
  internal::DelegatingArrayInputStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |queue_|.
  std::mutex lock_;
  std::condition_variable queue_has_room_;
  std::condition_variable queue_has_elements_;
  // The |queue_| contains the |Bytes| object filled by |Push| and not yet
  // consumed by |Pull|.  The two queues are out of step: an element is removed
  // from |queue_| by |Pull| when it returns a chunk to the stream, but the
  // corresponding element is removed from |done_| when |Pull| returns.  It is
  // executed at this point to destroy the chunk.
  std::queue<Bytes> queue_ GUARDED_BY(lock_);
  std::queue<std::function<void()>> done_ GUARDED_BY(lock_);
};

}  // namespace base
}  // namespace principia

#include "base/push_deserializer_body.hpp"
