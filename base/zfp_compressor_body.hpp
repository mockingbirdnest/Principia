#pragma once

#include "base/zfp_compressor.hpp"

namespace principia {
namespace base {
namespace zfp_compressor_internal {

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

  // Round up the size of the vector to a multiple of the block size.  This will
  // lead to poor compression at the end, but there is no support for "ignored"
  // data in zfp at this point.
  v.resize(((v.size() + block_ - 1) / block_) * block_, 0);
  CHECK_EQ(0, v.size() % block_);

  // Beware nx and ny!  (And the Jabberwock, my son!)
  // See https://zfp.readthedocs.io/en/release0.5.5/tutorial.html#high-level-c-interface
  std::unique_ptr<zfp_field, std::function<void(zfp_field*)>> const field(
      zfp_field_2d(v.data(),
                   /*type=*/zfp_type_double,
                   /*nx=*/block_,
                   /*ny=*/v.size() / block_),
      [](zfp_field* const field) { zfp_field_free(field); });
  WriteToMessage(field.get(), message);
}

template<int N>
void ZfpCompressor::ReadFromMessageNDimensional(
    std::vector<double>& v,
    std::string_view& message) const {
  if (v.empty()) {
    return;
  }

  // Make sure that we have enough space in the vector to decompress the
  // padding.
  v.resize(((v.size() + block_ - 1) / block_) * block_, 0);
  CHECK_EQ(0, v.size() % block_);

  std::unique_ptr<zfp_field, std::function<void(zfp_field*)>> const field(
      zfp_field_2d(v.data(),
                   /*type=*/zfp_type_double,
                   /*nx=*/block_,
                   /*ny=*/v.size() / block_),
      [](zfp_field* const field) { zfp_field_free(field); });

  ReadFromMessage(field.get(), message);
}

}  // namespace zfp_compressor_internal

using zfp_compressor_internal::ZfpCompressor;

}  // namespace base
}  // namespace principia
