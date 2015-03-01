#pragma once

#include <cstdint>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <memory>
#include <thread>

#include "base/not_null.hpp"

namespace principia {
namespace serialization {

class SynchronizingArrayOutputString
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  SynchronizingArrayOutputString(std::uint8_t* data, int const size);

  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;

 private:
  std::uint8_t* const data_;  // The byte array.
  const int size_;            // Total size of the array.

  int position_;
  int last_returned_size_;   // How many bytes we returned last time Next()
                             // was called (used for error checking only).
};

class Serializer {
 public:
  explicit Serializer(int const chunk_size);

  void Start(base::not_null<google::protobuf::Message const*> const message);

  base::not_null<std::string const*> Get();

private:
  std::unique_ptr<std::uint8_t[]> data_;
  SynchronizingArrayOutputString stream_;
  std::thread thread_;
};

}  // namespace serialization
}  // namespace principia