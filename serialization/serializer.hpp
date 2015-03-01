#pragma once

#include <cstdint>
#include <google/protobuf/message.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <thread>

#include "base/not_null.hpp"

namespace principia {
namespace serialization {

class SynchronizingArrayOutputString
    : public google::protobuf::io::ZeroCopyOutputStream {
 public:
  SynchronizingArrayOutputString(void* data, int size, int block_size = -1);
  bool Next(void** data, int* size) override;
  void BackUp(int count) override;
  std::int64_t ByteCount() const override;
};

class Serializer {
 public:
  explicit Serializer(int const chunk_size);

  void Start(base::not_null<google::protobuf::Message const*> const message);

  base::not_null<std::string const*> Get();

private:
  base::not_null<google::protobuf::Message const*> const message_;
  std::thread thread_;
};

}  // namespace serialization
}  // namespace principia