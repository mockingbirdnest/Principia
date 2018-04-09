
#pragma once

#include <cstdint>

#include "base/array.hpp"
#include "gipfeli/sinksource.h"

namespace principia {
namespace base {
namespace internal_sink_source {

template<typename Element>
class ArraySource : public google::compression::Source {
 public:
  explicit ArraySource(Array<Element> const& array);
  ~ArraySource() = default override;

  virtual size_t Available() const override;

  virtual const char* Peek(std::int64_t* length) override;

  virtual void Skip(std::int64_t n) override;

 private:
  const Array<Element> array_;
  std::int64_t next_to_read_ = 0;
};

template<typename Element>
class ArraySink : public google::compression::Sink {
 public:
  explicit ArraySink(Array<Element> const& array);
  ArraySink() = default override;

  virtual void Append(const char* data, size_t n) override;

  char* GetAppendBufferVariable(size_t min_size, size_t desired_size_hint,
                                char* scratch, size_t scratch_size,
                                size_t* allocated_size) override;

 private:
  const Array<Element> array_;
  std::int64_t next_to_write_ = 0;
};

}  // namespace internal_sink_source

using internal_sink_source::ArraySink;
using internal_sink_source::ArraySource;

}  // namespace base
}  // namespace principia

#include "base/sink_source_body.hpp"
