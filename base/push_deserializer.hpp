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

class MyStream : public google::protobuf::io::ZeroCopyInputStream {
 public:
  explicit MyStream(std::function<Bytes()> on_empty);
  ~MyStream() = default;

  bool Next(const void** data, int* size) override;
  void BackUp(int count) override;
  bool Skip(int count) override;
  std::int64_t ByteCount() const override;

 private:
  int size_;
  not_null<std::uint8_t const*> data_;
  std::function<Bytes()> on_empty_;

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

}  // namespace internal

class PushDeserializer {
 public:
  PushDeserializer(int const chunk_size,
                   int const number_of_chunks);
  ~PushDeserializer();

  void Start(not_null<google::protobuf::Message*> const message);

  void Push(Bytes const bytes);

 private:
  Bytes Pull();

  int const chunk_size_;
  int const number_of_chunks_;
  internal::MyStream stream_;
  std::unique_ptr<std::thread> thread_;

  // Synchronization objects for the |queue_|, which contains the |Bytes|
  // object filled by |Push| and not yet consumed by |Pull|.
  std::mutex lock_;
  std::condition_variable queue_has_room_;
  std::condition_variable queue_has_elements_;
  std::queue<Bytes> queue_ GUARDED_BY(lock_);
};

}  // namespace base
}  // namespace principia

#include "base/push_deserializer_body.hpp"
