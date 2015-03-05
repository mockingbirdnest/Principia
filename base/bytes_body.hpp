#pragma once

#include "base/bytes.hpp"

namespace principia {
namespace base {

Bytes::Bytes()
    : data(&null_data_), size(0) {}

Bytes::Bytes(base::not_null<std::uint8_t const*> const data,
             int const size)
    : data(data), size(size) {}

std::uint8_t Bytes::null_data_ = 0;

}  // namespace base
}  // namespace principia
