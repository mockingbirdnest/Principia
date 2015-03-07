#pragma once

#include "base/bytes.hpp"

namespace principia {
namespace base {

inline Bytes::Bytes()
    : data(reinterpret_cast<std::uint8_t*>(0xDEADBEEF)), size(0) {}

inline Bytes::Bytes(base::not_null<std::uint8_t*> const data,
                    int const size)
    : data(data), size(size) {}

}  // namespace base
}  // namespace principia
