#pragma once

#include <cstddef>
#include <cstdint>

#include "base/array.hpp"
#include "gipfeli/sinksource.h"

namespace principia {
namespace base {
namespace _sink_source {
namespace internal {

using namespace principia::base::_array;

template<typename Element>
class ArraySource : public google::compression::Source {
 public:
  explicit ArraySource(Array<Element> const& array);
  ~ArraySource() override = default;

  std::size_t Available() const override;

  const char* Peek(std::size_t* length) override;

  void Skip(std::size_t n) override;

 private:
  const Array<Element> array_;
  std::int64_t next_to_read_ = 0;
};

template<typename Element>
class ArraySink : public google::compression::Sink {
 public:
  explicit ArraySink(Array<Element> const& array);
  ~ArraySink() override = default;

  Array<Element> array() const;

  void Append(const char* data, std::size_t n) override;

  char* GetAppendBuffer(std::size_t min_size,
                        std::size_t desired_size_hint,
                        char* scratch,
                        std::size_t scratch_size,
                        std::size_t* allocated_size) override;

 private:
  const Array<Element> array_;
  std::int64_t next_to_write_ = 0;
};

}  // namespace internal

using internal::ArraySink;
using internal::ArraySource;

}  // namespace _sink_source
}  // namespace base
}  // namespace principia

#include "base/sink_source_body.hpp"
