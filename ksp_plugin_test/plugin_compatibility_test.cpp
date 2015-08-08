#include <fstream>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/coded_stream.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::Array;
using base::UniqueBytes;
using base::HexadecimalDecode;
using geometry::Bivector;
using geometry::Trivector;
using geometry::Vector;
using quantities::Length;
using si::Metre;

namespace ksp_plugin {

class PluginCompatibilityTest : public testing::Test {
};

TEST_F(PluginCompatibilityTest, PreBorel) {
  serialization::Multivector message;

  Vector<Length, Barycentric> const v({ -1 * Metre, 2 * Metre, 3 * Metre });
  v.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Vector<Length, Barycentric> const w =
    Vector<Length, Barycentric>::ReadFromMessage(message);
  Vector<Length, Barycentric> const expected_w(
  { -1 * Metre, 3 * Metre, 2 * Metre });
  EXPECT_EQ(expected_w, w);

  Bivector<Length, Barycentric> const b({ 4 * Metre, 5 * Metre, -6 * Metre });
  b.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Bivector<Length, Barycentric> const c =
    Bivector<Length, Barycentric>::ReadFromMessage(message);
  Bivector<Length, Barycentric> const expected_c(
  { -4 * Metre, 6 * Metre, -5 * Metre });
  EXPECT_EQ(expected_c, c);

  Trivector<Length, Barycentric> const t(-7 * Metre);
  t.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Trivector<Length, Barycentric> const u =
    Trivector<Length, Barycentric>::ReadFromMessage(message);
  Trivector<Length, Barycentric> const expected_u(7 * Metre);
  EXPECT_EQ(expected_u, u);
}

TEST_F(PluginCompatibilityTest, PreBourbaki) {
  std::fstream file = std::fstream("pre_bourbaki.proto.hex");
  CHECK(file.good());

  // Read the entire hex data.
  std::string hex;
  while (!file.eof()) {
    std::string line;
    std::getline(file, line);
    for (auto const c : line) {
      if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
        hex.append(1, c);
      }
    }
  }

  // Parse it and convert to binary data.
  UniqueBytes bin(hex.size() / 2);
  HexadecimalDecode(Array<std::uint8_t const>(
                        reinterpret_cast<const std::uint8_t*>(hex.c_str()),
                        hex.size()),
                    bin.get());

  // Construct a protocol buffer from the binary data.
  google::protobuf::io::CodedInputStream coded_input_stream(
      bin.data.get(), static_cast<int>(bin.size));
  serialization::Plugin serialized_plugin;
  CHECK(serialized_plugin.MergeFromCodedStream(&coded_input_stream));

  auto plugin = Plugin::ReadFromMessage(serialized_plugin);

  file.close();
}

}  // namespace ksp_plugin
}  // namespace principia