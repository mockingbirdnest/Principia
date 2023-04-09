#include <cstdint>
#include <string>

#include "astronomy/solar_system_fingerprints.hpp"
#include "astronomy/stabilize_ksp.hpp"
#include "base/serialization.hpp"
#include "geometry/frame.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace astronomy {

using ::testing::Eq;
using namespace principia::astronomy::_solar_system_fingerprints;
using namespace principia::astronomy::_stabilize_ksp;
using namespace principia::geometry::_frame;
using namespace principia::physics::_solar_system;

class KSPFingerprintTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<serialization::Frame::PluginTag,
                            Inertial,
                            Handedness::Right,
                            serialization::Frame::BARYCENTRIC>;

  KSPFingerprintTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt") {
    google::LogToStderr();
  }

  SolarSystem<Barycentric> solar_system_;
};

TEST_F(KSPFingerprintTest, Stock) {
  std::uint64_t const fingerprint = solar_system_.Fingerprint();
  LOG(INFO) << "Stock KSP fingerprint is 0x" << std::hex << std::uppercase
            << fingerprint;
  EXPECT_THAT(fingerprint, Eq(KSPStockSystemFingerprints[KSP191]));
}

TEST_F(KSPFingerprintTest, Corrected) {
  StabilizeKSP(solar_system_);
  std::uint64_t const fingerprint = solar_system_.Fingerprint();
  LOG(INFO) << "Corrected KSP fingerprint is 0x" << std::hex << std::uppercase
            << fingerprint;
  EXPECT_THAT(fingerprint, Eq(KSPStabilizedSystemFingerprints[KSP191]));
}

}  // namespace astronomy
}  // namespace principia
