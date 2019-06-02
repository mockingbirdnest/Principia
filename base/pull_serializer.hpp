
#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <queue>
#include <thread>

#include "absl/base/thread_annotations.h"
#include "absl/synchronization/mutex.h"
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
  DelegatingArrayOutputStream(
      Array<std::uint8_t> bytes,
      std::function<Array<std::uint8_t>(Array<std::uint8_t> bytes)> on_full);

  // The ZeroCopyOutputStream API.
  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  Array<std::uint8_t> bytes_;
  std::function<Array<std::uint8_t>(Array<std::uint8_t> bytes)> on_full_;

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
  // The |size| of the data objects enqueued by |Push| is never greater than
  // |chunk_size|.  At most |number_of_chunks| chunks are held in the internal
  // queue.  This class uses at most
  // |number_of_chunks * (chunk_size + O(1)) + O(1)| bytes.  Note that in the
  // presence of compression |chunk_size| is replaced by |compressed_chunk_size|
  // in this formula.
  PullSerializer(int chunk_size,
                 int number_of_chunks,
                 std::unique_ptr<Compressor> compressor);
  ~PullSerializer();

  // Starts the serializer, which will proceed to serialize |message|.  This
  // method must be called at most once for each serializer object.
  void Start(
      not_null<std::unique_ptr<google::protobuf::Message const>> message);
  void Start(not_null<google::protobuf::Message const*> message);

  // Obtain the next chunk of data from the serializer.  Blocks if no data is
  // available.  Returns a |Array<std::uint8_t>| object of |size| 0 at the end
  // of the serialization.  The returned object may become invalid the next time
  // |Pull| is called.
  // In the absence of compression, the data produced by |Pull| constitute a
  // stream and the boundaries between chunks are irrelevant.  In the presence
  // of compression however, the data producted by |Pull| are made of blocks and
  // the boundaries between chunks are relevant and must be preserved by the
  // clients and used when feeding data back to the deserializer.
  Array<std::uint8_t> Pull();

 private:
  // Enqueues the chunk of data to be returned to |Pull| and returns a free
  // chunk.  Blocks if there are no free chunks.  Used as a callback for the
  // underlying |DelegatingArrayOutputStream|.
  Array<std::uint8_t> Push(Array<std::uint8_t> bytes);

  // |owned_message_| is null if this object doesn't own the message.
  // |message_| is non-null after Start.
  std::unique_ptr<google::protobuf::Message const> owned_message_;
  google::protobuf::Message const* message_ = nullptr;

  std::unique_ptr<Compressor> const compressor_;

  // The chunk size passed at construction.  The stream outputs chunks of that
  // size.
  int const chunk_size_;

  // The maximum size of a chunk after compression.  Greater than |chunk_size_|
  // because the compressor will occasionally expand data.  This is the size of
  // the chunks in |data_|.
  int const compressed_chunk_size_;

  // The number of chunks passed at construction, used to size |data_|.
  int const number_of_chunks_;

  // How many of the |number_of_chunks_| chunks in |data_| are reserved for
  // compression.
  int const number_of_compression_chunks_;

  // The array supporting the stream and the stream itself.
  std::unique_ptr<std::uint8_t[]> data_;
  DelegatingArrayOutputStream stream_;

  // The thread doing the actual serialization.
  std::unique_ptr<std::thread> thread_;

  absl::Mutex lock_;

  // The |queue_| contains the |Array<std::uint8_t>| objects filled by |Push|
  // and not yet consumed by |Pull|.  If a |Array<std::uint8_t>| object has been
  // handed over to the caller by |Pull| it stays in the queue until the next
  // call to |Pull|, to make sure that the pointer is not reused while the
  // caller processes it.
  std::queue<Array<std::uint8_t>> queue_ GUARDED_BY(lock_);

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
