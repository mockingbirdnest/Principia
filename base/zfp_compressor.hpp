#pragma once

#include <optional>
#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "zfp/zfp.h"

namespace principia {
namespace base {
namespace zfp_compressor_internal {

using base::not_null;

// Helper class for ZFP compression.  This class compresses doubles, not
// quantities, because we don't want to depend on the layout of Quantity in
// memory.  Removal of unit dimensions must be done by the caller.
// Note that this is not a google::compression::Compressor, but it could be made
// into one.
class ZfpCompressor {
 public:
  // A compressor created with this constructor can only be used to read, not to
  // write.
  ZfpCompressor() = default;

  // A compressor created with this constructor can both read and write.
  explicit ZfpCompressor(double accuracy);

  // Read/write the version of ZFP from/to the message.
  template<typename Message>
  static void WriteVersion(not_null<Message*> message);
  template<typename Message>
  static void ReadVersion(Message const& message);

  // Serialization/deserialization of a vector of double into a message using an
  // N-dimensional encoding.  This encodes blocks of 4^N doubles and therefore
  // takes advantage of any correlation across these doubles.  The vector |v|
  // may be modified by padding it with zeroes.  When reading, the |message|
  // parameter is updated to reflect the data that was consumed.
  template<int N>
  void WriteToMessageNDimensional(std::vector<double>& v,
                                  not_null<std::string*> message) const;
  template<int N>
  void ReadFromMessageNDimensional(std::vector<double>& v,
                                   std::string_view& message) const;

  // Low-level API: serialization/deserialization of a field (allocated and
  // owned by the caller) into a message (which is expected to by a bytes field
  // of a proto).  When reading, the |message| parameter is updated to reflect
  // the data that was consumed.
  void WriteToMessage(const zfp_field* field,
                      not_null<std::string*> message) const;
  void ReadFromMessage(zfp_field* field,
                       std::string_view& message) const;

 private:
  static constexpr int block_ = 4;

  std::optional<double> const accuracy_;
};

}  // namespace zfp_compressor_internal

using zfp_compressor_internal::ZfpCompressor;

}  // namespace base
}  // namespace principia

#include "base/zfp_compressor_body.hpp"
