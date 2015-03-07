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

//TODO(phl): Revise all the comments.
namespace principia {
namespace base {

namespace internal {
// An output stream based on two arrays that delegates to a function the
// handling of the case where one array is full.  It calls the |on_full|
// function passed at construction and proceeds with filling the other array.
class DelegatingArrayOutputStream
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  // The stream is supported by |data| which has size |size << 1|.  Internally,
  // |data| is split into two pieces of size |size| and |on_full| is called when
  // one of the pieces is full.  Output can continue with the other piece.
  // |on_full| is expected to somehow consume the data of the array.  Note that
  // two successive calls to |on_full| are always passed different pointers, so
  // the caller can assume that it's fine to operate on that pointer in between
  // two calls to |on_full|.
  DelegatingArrayOutputStream(Bytes const bytes,
                              std::function<Bytes(Bytes const bytes)> on_full);

  // The ZeroCopyOutputStream API.
  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  int size_;
  not_null<std::uint8_t*> data_;
  std::function<Bytes(Bytes const bytes)> on_full_;

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

}  // namespace internal

// This class support serialization which is "pulled" by the client.  That is,
// the client creates a |PullSerializer| object, calls |Start| to start the
// serialization process, repeatedly calls |Pull| to obtain a chunk of data, and
// finally destroys the |PullSerializer|.  |PullSerializer| is intended for use
// in memory-critical contexts as it bounds the amount of memory used
// irrespective of the size of the message to serialize.
class PullSerializer {
 public:
  // The |size| of the data objects returned by |Pull| are never greater than
  // |max_size|.  This class uses at most |2 * max_size| bytes (plus some small
  // overhead).
  PullSerializer(int const chunk_size, int const number_of_chunks);
  ~PullSerializer();

  // Starts the serializer, which will proceed to serialize |message|.  This
  // method must be called at most once for each serializer object.
  void Start(not_null<google::protobuf::Message const*> const message);

  // Obtain the next chunk of data from the serializer.  Blocks if no data is
  // available.  Returns a |Bytes| object of |size| 0 at the end of the
  // serialization.
  Bytes Pull();

 private:
  // Sets the chunk of data to be returned to |Pull|.  Used as a callback for
  // the underlying |DelegatingTwoArrayOutputStream|.
  Bytes Push(Bytes const bytes);

  int const chunk_size_;
  int const number_of_chunks_;

  std::unique_ptr<std::uint8_t[]> data_;
  internal::DelegatingArrayOutputStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |queue_|, which contains the |Bytes|
  // object filled by |Push| and not yet consumed by |Pull|.
  std::mutex lock_;
  std::condition_variable queue_has_room_;
  std::condition_variable queue_has_elements_;
  std::queue<Bytes> queue_ GUARDED_BY(lock_);
  std::queue<not_null<std::uint8_t*>> free_ GUARDED_BY(lock_);
  bool is_first_pull_ GUARDED_BY(lock_);
};

}  // namespace base
}  // namespace principia

#include "base/pull_serializer_body.hpp"
