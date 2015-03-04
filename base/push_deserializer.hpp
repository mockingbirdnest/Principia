#pragma once

#include <cstdint>
#include <condition_variable>  // NOLINT(build/c++11)
#include <memory>
#include <mutex>  // NOLINT(build/c++11)
#include <thread>  // NOLINT(build/c++11)

#include "base/bytes.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "google/protobuf/message.h"
#include "google/protobuf/io/zero_copy_stream.h"

namespace principia {
namespace base {

namespace internal {

class MyStream : public google::protobuf::io::ZeroCopyInputStream {
 public:
  MyStream(not_null<std::uint8_t*> data,
           int const size,
           std::function<Bytes()> on_empty);
  ~MyStream() = default;

  bool Next(const void** data, int* size) override;
  void BackUp(int count) override;
  bool Skip(int count) override;
  std::int64_t ByteCount() const override;

 private:
  int const size_;
  not_null<std::uint8_t*> data1_;
  not_null<std::uint8_t*> data2_;
  std::function<Bytes()> on_empty_;

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

}  // namespace internal

class PushDeserializer {
 public:
  explicit PushDeserializer(int const max_size);
  ~PushDeserializer();

  void Start(not_null<google::protobuf::Message*> const message);

  void Push(Bytes const bytes);

 private:
  Bytes Pull();

  std::unique_ptr<std::uint8_t[]> data_;
  internal::MyStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |holder_|, which contains the |Bytes|
  // object filled by |Push| and not yet consumed by |Pull|.  The |holder_| is
  // effectively a 1-element queue.
  std::mutex lock_;
  std::condition_variable holder_is_empty_;
  std::condition_variable holder_is_full_;
  std::unique_ptr<Bytes> holder_ GUARDED_BY(lock_);  // std::optional, really.

  // Set to true when the |thread_| completes.
  bool done_ GUARDED_BY(lock_);
};

}  // namespace base
}  // namespace principia

#include "base/push_deserializer_body.hpp"
