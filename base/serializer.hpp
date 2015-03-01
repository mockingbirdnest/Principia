#pragma once

#include <cstdint>
#include <condition_variable>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <memory>
#include <mutex>
#include <thread>

#include "base/macros.hpp"
#include "base/not_null.hpp"

namespace principia {
namespace base {

// An output stream based on two arrays that delegates to a function the
// handling of the case where one array is full.  It calls the |on_full|
// function passed at construction and proceeds with filling the other array.
class DelegatingTwoArrayOutputStream
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  // The stream is supported by |data| which has size |size|.  Internally,
  // |data| is split into two pieces and |on_full| is called when one of the
  // pieces is full.  Output can continue with the other piece.  |on_full| is
  // expected to somehow consume the data of the array.  Note that two
  // successive calls to |on_full| are always passed different pointers, so the
  // caller can assume that it's fine to operate on that pointer in between two
  // calls to |on_full|.
  DelegatingTwoArrayOutputStream(
      not_null<std::uint8_t*> data,
      int const size,
      std::function<void(not_null<std::uint8_t const*> const data,
                         int const size)> on_full);

  // The ZeroCopyOutputStream API.
  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  int const size_;
  not_null<std::uint8_t*> data1_;
  not_null<std::uint8_t*> data2_;
  std::function<void(not_null<std::uint8_t const*> const data,
                     int const size)> on_full_;

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

class PullSerializer {
 public:
  // Similar to a StringPiece.  |data| is not owned.
  struct Data {
    Data(not_null<std::uint8_t const*> const data, int const size);
    not_null<std::uint8_t const*> const data;
    int const size;
  };

  // The |size| of the data objects returned by |Get| are never greater than
  // |max_size|.
  explicit PullSerializer(int const max_size);
  ~PullSerializer();

  // Starts the serializer, which will proceed to serialize |message|.  This
  // method must be called at most once for each serializer object.
  void Start(not_null<google::protobuf::Message const*> const message);

  // Obtain the next chunk of data from the serializer.  Blocks if no data is
  // available.  Returns a |Data| object of |size| 0 at the end of the
  // serialization.
  Data Get();

private:
  // Sets the chunk of data to be returned to |Get|.  Used as a callback for the
  // underlying |DelegatingTwoArrayOutputStream|.
  void Set(not_null<std::uint8_t const*> const data, int const size);

  // Placeholders for an empty |Data| object.
  static std::uint8_t no_data_data_;
  static Data* no_data_;

  std::unique_ptr<std::uint8_t[]> data_;
  DelegatingTwoArrayOutputStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |holder_|, which contains the |Data|
  // object filled by |Set| and not yet consumed by |Get|.
  std::mutex lock_;
  std::condition_variable holder_is_empty_;
  std::condition_variable holder_is_full_;
  std::unique_ptr<Data> holder_ GUARDED_BY(lock_);  // std::optional, really.

  // Set to true when the |thread_| completes.
  bool done_ GUARDED_BY(lock_);
};

}  // namespace base
}  // namespace principia

#include "base/PullSerializer_body.hpp"
