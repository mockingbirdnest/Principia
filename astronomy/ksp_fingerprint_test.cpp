
#include <string>

#include "astronomy/solar_system_fingerprints.hpp"
#include "astronomy/stabilize_ksp.hpp"
#include "base/fingerprint2011.hpp"
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

using base::Fingerprint2011;
using base::SerializeAsBytes;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using physics::SolarSystem;
using ::testing::Eq;

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
  auto const hierarchical_system = solar_system_.MakeHierarchicalSystem();
  serialization::HierarchicalSystem message;
  hierarchical_system->WriteToMessage(&message);
  uint64_t const fingerprint = Fingerprint2011(SerializeAsBytes(message).get());
  LOG(INFO) << "Stock KSP fingerprint is 0x" << std::hex << std::uppercase
            << fingerprint;
  EXPECT_THAT(fingerprint, Eq(KSPStockSystemFingerprints[KSP191]));
}

TEST_F(KSPFingerprintTest, Corrected) {
  StabilizeKSP(solar_system_);
  auto const hierarchical_system = solar_system_.MakeHierarchicalSystem();
  serialization::HierarchicalSystem message;
  hierarchical_system->WriteToMessage(&message);
  uint64_t const fingerprint = Fingerprint2011(SerializeAsBytes(message).get());
  LOG(INFO) << "Corrected KSP fingerprint is 0x" << std::hex << std::uppercase
            << fingerprint;
  EXPECT_THAT(fingerprint, Eq(KSPStabilizedSystemFingerprint[KSP191]));
}

}  // namespace astronomy
}  // namespace principia
