#pragma once

#include <cstdint>
#include <condition_variable>  // NOLINT(build/c++11)
#include <memory>
#include <mutex>  // NOLINT(build/c++11)
#include <thread>  // NOLINT(build/c++11)

#include "base/bytes.hpp"
#include "base/not_null.hpp"
#include "google/protobuf/message.h"
#include "google/protobuf/io/zero_copy_stream.h"

namespace principia {
namespace base {

namespace internal {

class MyStream : public google::protobuf::io::ZeroCopyInputStream {
 public:
  MyStream(const void* data, int size, int block_size = -1);
  ~MyStream() = default;

  bool Next(const void** data, int* size) override;
  void BackUp(int count) override;
  bool Skip(int count) override;
  std::int64_t ByteCount() const override;

 private:
  const std::uint8_t* const data_;  // The byte array.
  const int size_;           // Total size of the array.
  const int block_size_;     // How many bytes to return at a time.

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

}  // namespace internal

class PushDeserializer {
 public:
  explicit PushDeserializer(int const max_size);
  ~PushDeserializer();

 private:
  std::unique_ptr<std::uint8_t[]> data_;
  internal::MyStream stream_;
  std::unique_ptr<std::thread> thread_;

};

}  // namespace base
}  // namespace principia

#include "base/push_deserializer_body.hpp"
