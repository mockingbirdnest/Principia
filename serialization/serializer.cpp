#include "serialization/serializer.hpp"

using std::placeholders::_1;
using std::placeholders::_2;
using std::swap;

namespace principia {
namespace serialization {

namespace {

void Serialize(base::not_null<google::protobuf::Message const*> const message,
               google::protobuf::io::ZeroCopyOutputStream* const stream) {
  CHECK(message->SerializeToZeroCopyStream(stream));
}

}  // namespace

SynchronizingArrayOutputString::SynchronizingArrayOutputString(
    base::not_null<std::uint8_t*> data,
    int const size,
    std::function<void(base::not_null<std::uint8_t const*> const data,
                       int const size)> on_full)
    : size_(size >> 1),
      data1_(&data[0]),
      data2_(&data[size_]),
      on_full_(std::move(on_full)),
      position_(0),
      last_returned_size_(0) {
  CHECK_EQ(0, size % 2);
}

bool SynchronizingArrayOutputString::Next(void** data, int* size) {
  if (position_ < size_) {
    last_returned_size_ = size_ - position_;
    *data = &data1_[position_];
    *size = last_returned_size_;
    position_ += last_returned_size_;
  } else {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the other array.
    on_full_(data1_, size_);
    position_ = 0;
    last_returned_size_ = 0;
    swap(data1_, data2_);
  }
  return true;
}

void SynchronizingArrayOutputString::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  position_ -= count;
  last_returned_size_ = 0;  // Don't let caller back up further.
}

std::int64_t SynchronizingArrayOutputString::ByteCount() const {
  return position_;
}

Serializer::Data::Data(base::not_null<std::uint8_t const*> const data,
                       int const size)
    : data(data), size(size) {}

Serializer::Serializer(int const chunk_size)
    : data_(std::make_unique<std::uint8_t[]>(chunk_size)),
      stream_(data_.get(),
              chunk_size,
              std::bind(&Serializer::Set, this, _1, _2)) {}

Serializer::~Serializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

void Serializer::Start(
    base::not_null<google::protobuf::Message const*> const message) {
  thread_ = std::make_unique<std::thread>([this, message](){
    CHECK(message->SerializeToZeroCopyStream(&stream_));
  });
}

Serializer::Data Serializer::Get() {
  Data* result;
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_full_.wait(l, [this](){ return holder_ != nullptr; });
    result = holder_.release();
  }
  holder_is_empty_.notify_all();
  return *result;
}

void Serializer::Set(base::not_null<std::uint8_t const*> const data,
                     int const size) {
  {
    std::unique_lock<std::mutex> l(lock_);
    holder_is_empty_.wait(l, [this](){ return holder_ == nullptr; });
    holder_ = std::make_unique<Data>(data, size);
  }
  holder_is_full_.notify_all();
}

}  // namespace serialization
}  // namespace principia
