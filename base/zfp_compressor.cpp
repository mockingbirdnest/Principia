#include "base/zfp_compressor.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "base/array.hpp"

namespace principia {
namespace base {
namespace zfp_compressor_internal {

ZfpCompressor::ZfpCompressor(double const accuracy) : accuracy_(accuracy) {}

void ZfpCompressor::WriteToMessage2D(
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

void ZfpCompressor::ReadFromMessage2D(std::vector<double>& v,
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

  zfp_write_header(zfp.get(), field, ZFP_HEADER_FULL);
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
  size_t const header_bits = zfp_read_header(zfp.get(), field, ZFP_HEADER_FULL);
  CHECK_LT(0, header_bits);

  size_t const compressed_size = zfp_decompress(zfp.get(), field);
  CHECK_LT(0, compressed_size);
  message.remove_prefix(compressed_size);
}

}  // namespace zfp_compressor_internal
}  // namespace base
}  // namespace principia
