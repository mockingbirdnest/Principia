﻿
#include <filesystem>
#include <fstream>
#include <map>
#include <string>

#include "base/allocated_by.hpp"
#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/coded_stream.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using base::Array;
using base::UniqueArray;
using base::AssertAllocatedBy;
using base::HexadecimalEncoder;
using geometry::Bivector;
using geometry::Trivector;
using geometry::Vector;
using quantities::Length;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::Gt;
using ::testing::Lt;

class TestablePlugin : public Plugin {
 public:
  void KeepAllVessels();

  std::map<GUID, not_null<Vessel const*>> vessels() const;

  static not_null<std::unique_ptr<TestablePlugin>> ReadFromMessage(
    serialization::Plugin const& message);
};

void TestablePlugin::KeepAllVessels() {
  for (auto const& [_, vessel] : vessels_) {
    kept_vessels_.insert(vessel.get());
  }
}

std::map<GUID, not_null<Vessel const*>> TestablePlugin::vessels() const {
  std::map<GUID, not_null<Vessel const*>> result;
  for (auto const& [guid, vessel] : vessels_) {
    result.insert(std::make_pair(guid, vessel.get()));
  }
  return result;
}

not_null<std::unique_ptr<TestablePlugin>> TestablePlugin::ReadFromMessage(
  serialization::Plugin const& message) {
  Plugin* plugin = Plugin::ReadFromMessage(message).release();
  return std::unique_ptr<TestablePlugin>(
      AssertAllocatedBy<std::allocator<TestablePlugin>>(
          static_cast<TestablePlugin*>(plugin)));
}

class PluginCompatibilityTest : public testing::Test {
 protected:
  serialization::Plugin ReadFromFile(std::string const& filename) {
    // Open the file and read hexadecimal data.
    std::fstream file =
        std::fstream(SOLUTION_DIR / "ksp_plugin_test" / filename);
    CHECK(file.good());
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
    file.close();

    // Parse the hexadecimal data and convert it to binary data.
    HexadecimalEncoder</*null_terminated=*/false> encoder;
    auto const bin = encoder.Decode(hex);

    // Construct a protocol buffer from the binary data.
    google::protobuf::io::CodedInputStream coded_input_stream(
        bin.data.get(), static_cast<int>(bin.size));
    serialization::Plugin message;
    CHECK(message.MergeFromCodedStream(&coded_input_stream));

    return message;
  }
};

TEST_F(PluginCompatibilityTest, PreCartan) {
  // This space for rent.
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
