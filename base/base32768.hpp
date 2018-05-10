
#pragma once

#include <cstdint>

#include "base/array.hpp"

namespace principia {
namespace base {
namespace internal_base32768 {

inline void Base32768Encode(Array<std::uint8_t const> input,
                            Array<char16_t> output);

// Same as above but the storage is allocated by the callee.  If
// |null_terminated| is true a null character is appended to the encoded form.
inline UniqueArray<char16_t> Base32768Encode(Array<std::uint8_t const> input,
                                             bool null_terminated);

inline void Base32768Decode(Array<char16_t const> input,
                            Array<std::uint8_t> output);

// Same as above but the storage is allocated by the callee.  The input may or
// may not be null-terminated.
inline UniqueArray<std::uint8_t> Base32768Decode(Array<char16_t const> input);

}  // namespace internal_base32768

using internal_base32768::Base32768Decode;
using internal_base32768::Base32768Encode;

}  // namespace base
}  // namespace principia

#include "base/base32768_body.hpp"
