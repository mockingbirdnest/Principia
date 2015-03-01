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
namespace serialization {

class SynchronizingArrayOutputString
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  SynchronizingArrayOutputString(
      base::not_null<std::uint8_t*> data,
      int const size,
      std::function<void(base::not_null<std::uint8_t const*> const data,
                         int const size)> on_full);

  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  int const size_;
  base::not_null<std::uint8_t*> data1_;
  base::not_null<std::uint8_t*> data2_;
  std::function<void(base::not_null<std::uint8_t const*> const data,
                     int const size)> on_full_;

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

class Serializer {
 public:
  struct Data {
    Data(base::not_null<std::uint8_t const*> const data, int const size);
    base::not_null<std::uint8_t const*> const data;
    int const size;
  };

  explicit Serializer(int const chunk_size);
  ~Serializer();

  void Start(base::not_null<google::protobuf::Message const*> const message);

  Data Get();

private:

  void Set(base::not_null<std::uint8_t const*> const data, int const size);

  std::unique_ptr<std::uint8_t[]> data_;
  SynchronizingArrayOutputString stream_;
  std::unique_ptr<std::thread> thread_;

  std::mutex lock_;
  std::condition_variable holder_is_empty_;
  std::condition_variable holder_is_full_;
  std::unique_ptr<Data> holder_ GUARDED_BY(lock_);
};

}  // namespace serialization
}  // namespace principia