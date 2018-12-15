
#pragma once

#include "base/base64.hpp"

#include <string>

#include "absl/strings/escaping.h"

namespace principia {
namespace base {
namespace internal_base64 {

constexpr std::int64_t bits_per_byte = 8;
constexpr std::int64_t bits_per_char = 6;

template<bool null_terminated>
void principia::base::internal_base64::Base64Encoder<null_terminated>::Encode(
    Array<std::uint8_t const> input,
    Array<char> output) {
  std::string_view const input_view(reinterpret_cast<const char*>(input.data),
                                    input.size);
  std::string output_string;
  absl::WebSafeBase64Escape(input_view, &output_string);
  std::memcpy(output.data, output_string.c_str(), output_string.size());
  if constexpr (null_terminated) {
    output.data[output_size.size()] = 0;
  }
}

template<bool null_terminated>
UniqueArray<char> Base64Encoder<null_terminated>::Encode(
    Array<std::uint8_t const> input) {
  UniqueArray<char> output(EncodedLength(input));
  if (output.size > 0) {
    Encode(input, output.get());
  }
  return output;
}

template<bool null_terminated>
std::int64_t Base64Encoder<null_terminated>::EncodedLength(
    Array<std::uint8_t const> input) {
  std::int64_t const nonterminated_length =
      (bits_per_byte * input.size + bits_per_char - 1) / bits_per_char;
  if constexpr (null_terminated) {
    return nonterminated_length + 1;
  } else {
    return nonterminated_length;
  }
}

template<bool null_terminated>
void Base64Encoder<null_terminated>::Decode(Array<char const> input,
                                            Array<std::uint8_t> output) {
  std::string_view const input_view(input.data, input.size);
  std::string output_string;
  absl::WebSafeBase64Unescape(input_view, &output_string);
  std::memcpy(output.data, output_string.c_str(), output_string.size());
}

template<bool null_terminated>
UniqueArray<std::uint8_t> Base64Encoder<null_terminated>::Decode(
    Array<char const> input) {
  UniqueArray<std::uint8_t> output(DecodedLength(input));
  if (output.size > 0) {
    Decode(input, output.get());
  }
  return output;
}

template<bool null_terminated>
std::int64_t Base64Encoder<null_terminated>::DecodedLength(
    Array<char const> input) {
  return bits_per_char * input.size / bits_per_byte;
}

}  // namespace internal_base64
}  // namespace base
}  // namespace principia
