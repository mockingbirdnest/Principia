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

}  // namespace zfp_compressor_internal

using zfp_compressor_internal::ZfpCompressor;

}  // namespace base
}  // namespace principia
