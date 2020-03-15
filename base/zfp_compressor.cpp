#include "base/zfp_compressor.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "base/array.hpp"

namespace principia {
namespace base {
namespace zfp_compressor_internal {

// ZFP headers limit the dimensions to 2^(48 / N).  For N = 4, this is way too
// small for our purposes.  Therefore, we must not include the bit
// ZFP_HEADER_META.  We still need a header, though, to record the compression
// parameters.
static constexpr uint header_mask = ZFP_HEADER_MODE | ZFP_HEADER_MAGIC;

ZfpCompressor::ZfpCompressor(double const accuracy) : accuracy_(accuracy) {}

void ZfpCompressor::WriteToMessage(const zfp_field* const field,
                                   not_null<std::string*> message) const {
  std::unique_ptr<zfp_stream, std::function<void(zfp_stream*)>> const zfp(
      zfp_stream_open(/*stream=*/nullptr),
      [](zfp_stream* const zfp) { zfp_stream_close(zfp); });

  CHECK(accuracy_.has_value());
  if (accuracy_ == 0) {
    zfp_stream_set_reversible(zfp.get());
  } else {
    zfp_stream_set_accuracy(zfp.get(), *accuracy_);
  }
  size_t const buffer_size = zfp_stream_maximum_size(zfp.get(), field);
  UniqueArray<std::uint8_t> const buffer(buffer_size);
  not_null<bitstream*> const stream =
      check_not_null(stream_open(buffer.data.get(), buffer_size));
  zfp_stream_set_bit_stream(zfp.get(), &*stream);

  zfp_write_header(zfp.get(), field, header_mask);
  size_t const compressed_size = zfp_compress(zfp.get(), field);
  CHECK_LT(0, compressed_size);
  message->append(static_cast<char const*>(stream_data(stream)),
                  stream_size(stream));
}

void ZfpCompressor::ReadFromMessage(zfp_field* const field,
                                    std::string_view& message) const {
  std::unique_ptr<zfp_stream, std::function<void(zfp_stream*)>> const zfp(
      zfp_stream_open(/*stream=*/nullptr),
      [](zfp_stream* const zfp) { zfp_stream_close(zfp); });

  not_null<bitstream*> const stream = check_not_null(
      stream_open(const_cast<char*>(&message.front()), message.size()));
  zfp_stream_set_bit_stream(zfp.get(), &*stream);
  size_t const header_bits = zfp_read_header(zfp.get(), field, header_mask);
  CHECK_LT(0, header_bits);

  size_t const compressed_size = zfp_decompress(zfp.get(), field);
  CHECK_LT(0, compressed_size);
  message.remove_prefix(compressed_size);
}

}  // namespace zfp_compressor_internal
}  // namespace base
}  // namespace principia
