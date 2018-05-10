
#pragma once

#include <cstdint>

#include "base/array.hpp"

// This file implements the base 32768 encoding defined by
// https://github.com/qntm/base32768 (MIT License).  This is a complete
// reimplementation in C++.

namespace principia {
namespace base {
namespace internal_base32768 {

// Encodes |input| into |output|, which must be large enough to hold the encoded
// form.
inline void Base32768Encode(Array<std::uint8_t const> input,
                            Array<char16_t> output);

// Same as above but the storage is allocated by the callee.  If
// |null_terminated| is true a null character is appended to the encoded form.
inline UniqueArray<char16_t> Base32768Encode(Array<std::uint8_t const> input,
                                             bool null_terminated);

// Length of the encoded form, in char16_t.
inline std::int64_t Base32768EncodedLength(Array<std::uint8_t const> input);

// Decodes |input| into |output|, which must be large enough to hold the decoded
// form.
inline void Base32768Decode(Array<char16_t const> input,
                            Array<std::uint8_t> output);

// Same as above but the storage is allocated by the callee.  The input may or
// may not be null-terminated.
inline UniqueArray<std::uint8_t> Base32768Decode(Array<char16_t const> input);

// Length of the decoded form, in uint8_t.
inline std::int64_t Base32768DecodedLength(Array<char16_t const> input);

}  // namespace internal_base32768

using internal_base32768::Base32768Decode;
using internal_base32768::Base32768DecodedLength;
using internal_base32768::Base32768Encode;
using internal_base32768::Base32768EncodedLength;

}  // namespace base
}  // namespace principia

#include "base/base32768_body.hpp"
