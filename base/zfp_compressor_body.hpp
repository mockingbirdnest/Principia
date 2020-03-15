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

  static unique_zfp_field NewField(zfp_type type, std::vector<double>& v);

 private:
  // Returns the size rounded up to a multiple of 4^(N - 1).
  static constexpr std::int64_t RoundUp(std::int64_t size);

  // The ZFP block size.
  static constexpr int block_ = 4;
  // Data must be padded to a multiple of that value, which is 4^(N - 1).
  static constexpr int padding_ = 1 << (2 * (N - 1));
};

template<int N>
typename NDimensionalHelper<N>::unique_zfp_field
NDimensionalHelper<N>::NewField(zfp_type const type, std::vector<double>& v) {
  auto free = [](zfp_field* const field) { zfp_field_free(field); };

  // On compression: Round up the size of the vector to a multiple of the block
  // size.  This will lead to poor compression at the end, but there is no
  // support for "ignored" data in zfp at this point.
  // On decompression: Make sure that we have enough space in the vector to
  // decompress the padding.
  v.resize(RoundUp(v.size()), 0);
  CHECK_EQ(0, v.size() % padding_);

  // Beware nx, ny and friends!  (And the Jabberwock, my son!)
  // See
  // https://zfp.readthedocs.io/en/release0.5.5/tutorial.html#high-level-c-interface
  if constexpr (N == 1) {
    return unique_zfp_field(
        zfp_field_1d(v.data(), type, /*nx=*/v.size() / padding_),
        std::move(free));
  } else if constexpr (N == 2) {
    return unique_zfp_field(
        zfp_field_2d(v.data(), type, /*nx=*/block_, /*ny=*/v.size() / padding_),
        std::move(free));
  } else if constexpr (N == 3) {
    return unique_zfp_field(zfp_field_3d(v.data(),
                                         type,
                                         /*nx=*/block_,
                                         /*ny=*/block_,
                                         /*nz=*/v.size() / padding_),
                            std::move(free));
  } else if constexpr (N == 4) {
    return unique_zfp_field(zfp_field_4d(v.data(),
                                         type,
                                         /*nx=*/block_,
                                         /*ny=*/block_,
                                         /*nz=*/block_,
                                         /*nw=*/v.size() / padding_),
                            std::move(free));
  } else {
    static_assert(false, "Unsupported dimension");
  }
}

template<int N>
constexpr std::int64_t NDimensionalHelper<N>::RoundUp(std::int64_t const size) {
  // Hacker's Delight, H. S. Warren, Jr., section 3-1.
  return (size + (padding_ - 1)) & (-padding_);
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

  auto const field =
      NDimensionalHelper<N>::NewField(/*type=*/zfp_type_double, v);
  WriteToMessage(field.get(), message);
}

template<int N>
void ZfpCompressor::ReadFromMessageNDimensional(
    std::vector<double>& v,
    std::string_view& message) const {
  if (v.empty()) {
    return;
  }

  auto const field =
      NDimensionalHelper<N>::NewField(/*type=*/zfp_type_double, v);
  ReadFromMessage(field.get(), message);
}

}  // namespace zfp_compressor_internal

using zfp_compressor_internal::ZfpCompressor;

}  // namespace base
}  // namespace principia
