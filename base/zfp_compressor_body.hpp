#pragma once

#include "base/zfp_compressor.hpp"

#include <functional>

namespace principia {
namespace base {
namespace zfp_compressor_internal {

template<int N>
class NDimensionalHelper {
 public:
  using unique_zfp_field =
      std::unique_ptr<zfp_field, std::function<void(zfp_field*)>>;

  // Returns the size rounded up to a multiple of 4^(N - 1).
  static constexpr std::int64_t RoundUp(std::int64_t size);

  static unique_zfp_field NewField(std::vector<double>& v,
                                   zfp_type type);

 private:
  // 4^(N - 1).
  static constexpr std::int64_t block_ = 1 << (2 * (N - 1));
};

template<int N>
constexpr std::int64_t NDimensionalHelper<N>::RoundUp(std::int64_t const size) {
  // Hacker's Delight, H. S. Warren, Jr., section 3-1.
  return (size + (block_ - 1)) & (-block_);
}

template<int N>
typename NDimensionalHelper<N>::unique_zfp_field
NDimensionalHelper<N>::NewField(std::vector<double>& v,
                                zfp_type const type) {
  // Beware nx and ny!  (And the Jabberwock, my son!)
  // See
  // https://zfp.readthedocs.io/en/release0.5.5/tutorial.html#high-level-c-interface
  CHECK_EQ(0, v.size() % block_);
  if constexpr (N == 1) {
  } else if constexpr (N == 2) {
    return unique_zfp_field(
        zfp_field_2d(v.data(), type, /*nx=*/block_, /*ny=*/v.size() / block_),
        [](zfp_field* const field) { zfp_field_free(field); });
  } else if constexpr (N == 3) {
  } else if constexpr (N == 4) {
  } else {
    static_assert(false, "Unsupported dimension");
  }
}

template<typename Message>
void ZfpCompressor::WriteVersion(not_null<Message*> message) {
  // For future compatibility, record the version of ZFP used to compress.
  auto* const zfp = message->mutable_zfp();
  zfp->set_codec_version(ZFP_CODEC);
  zfp->set_library_version(ZFP_VERSION);
}

template<typename Message>
void ZfpCompressor::ReadVersion(Message const& message) {
  // We only read for the fun of checking, for now.
  CHECK_EQ(ZFP_CODEC, message.zfp().codec_version());
  CHECK_EQ(ZFP_VERSION, message.zfp().library_version());
}

template<int N>
void ZfpCompressor::WriteToMessageNDimensional(
    std::vector<double>& v,
    not_null<std::string*> const message) const {
  if (v.empty()) {
    return;
  }
  using Helper = NDimensionalHelper<N>;

  // Round up the size of the vector to a multiple of the block size.  This will
  // lead to poor compression at the end, but there is no support for "ignored"
  // data in zfp at this point.
  v.resize(Helper::RoundUp(v.size()), 0);

  auto const field = Helper::NewField(v, /*type=*/zfp_type_double);
  WriteToMessage(field.get(), message);
}

template<int N>
void ZfpCompressor::ReadFromMessageNDimensional(
    std::vector<double>& v,
    std::string_view& message) const {
  if (v.empty()) {
    return;
  }
  using Helper = NDimensionalHelper<N>;

  // Make sure that we have enough space in the vector to decompress the
  // padding.
  v.resize(Helper::RoundUp(v.size()), 0);
  auto const field = Helper::NewField(v, /*type=*/zfp_type_double);
  ReadFromMessage(field.get(), message);
}

}  // namespace zfp_compressor_internal

using zfp_compressor_internal::ZfpCompressor;

}  // namespace base
}  // namespace principia
