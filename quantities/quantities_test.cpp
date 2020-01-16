
#include <functional>
#include <string>

#include "google/protobuf/stubs/common.h"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/astronomy.hpp"
#include "quantities/bipm.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "testing_utilities/algebra.hpp"

namespace principia {
namespace quantities {

using astronomy::JovianGravitationalParameter;
using astronomy::JulianYear;
using astronomy::LightYear;
using astronomy::Parsec;
using astronomy::TerrestrialGravitationalParameter;
using constants::ElectronMass;
using constants::ProtonMass;
using constants::SpeedOfLight;
using si::Ampere;
using si::Candela;
using si::Day;
using si::Hour;
using si::Joule;
using si::Kelvin;
using si::Kilogram;
using si::Metre;
using si::Mole;
using si::Radian;
using si::Second;
using si::Steradian;
using uk::Foot;
using uk::Furlong;
using uk::Gallon;
using uk::Mile;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::MatchesRegex;

class QuantitiesTest : public testing::Test {};

using QuantitiesDeathTest = QuantitiesTest;

TEST_F(QuantitiesTest, DimensionfulComparisons) {
  testing_utilities::TestOrder(TerrestrialGravitationalParameter,
                               JovianGravitationalParameter);
  testing_utilities::TestOrder(LightYear, Parsec);
  testing_utilities::TestOrder(-SpeedOfLight, SpeedOfLight);
  testing_utilities::TestOrder(SpeedOfLight * Day, LightYear);
}

TEST_F(QuantitiesTest, DimensionlfulOperations) {
  testing_utilities::TestVectorSpace(
      0 * Metre / Second, SpeedOfLight, 88 * Mile / Hour,
      -340.29 * Metre / Second, 0.0, 1.0, -2 * π, 1729.0, 0, 2);
  // Dimensionful multiplication is a tensor product, see [Tao 2012].
  testing_utilities::TestBilinearMap(
      std::multiplies<>(), ProtonMass, ElectronMass, SpeedOfLight,
      1 * Furlong / JulianYear, -e, 0, 2);
}

// The Greek letters cause a warning when stringified by the macros, because
// apparently Visual Studio doesn't encode strings in UTF-8 by default.
#pragma warning(disable: 4566)

TEST_F(QuantitiesTest, Formatting) {
  auto const all_the_units = 1 * Metre * Kilogram * Second * Ampere * Kelvin /
                                 (Mole * Candela * Radian * Steradian);
  std::string const expected = "+1e+00 m kg s A K mol^-1 cd^-1 rad^-3";
  std::string const actual = DebugString(all_the_units, 0);
  EXPECT_EQ(expected, actual);
  std::string const π17 = R"(\+3\.1415926535897931.e\+00)";
  EXPECT_THAT(DebugString(π), MatchesRegex(π17));
  std::string const minus_e17 = R"(\-2\.718281828459045..e\+00)";
  EXPECT_THAT(DebugString(-e), MatchesRegex(minus_e17));
}

#pragma warning(default: 4566)

TEST_F(QuantitiesTest, RotationalUnits) {
  EXPECT_THAT(SIUnit<AngularFrequency>(), Eq(Radian / Second));
  EXPECT_THAT(SIUnit<AngularAcceleration>(), Eq(Radian / Pow<2>(Second)));
  // SI Brochure 8th edition, 2006, updated in 2014, Section 2.2.2:
  // For example, the quantity torque may be thought of as the cross product of
  // force and distance, suggesting the unit newton metre, or it may be thought
  // of as energy per angle, suggesting the unit joule per radian.
  // But we do things differently, see the comment in the declaration of
  // MomentOfInertia.
  EXPECT_THAT(SIUnit<Torque>(), Eq(Joule * Radian));
}

TEST_F(QuantitiesTest, IsFinite) {
  Length l = 0 * Foot;
  // The serialization is to defeat the optimizer which tries to prevent us from
  // dividing by 0.
  serialization::Quantity message;
  l.WriteToMessage(&message);
  l = Length::ReadFromMessage(message);
  EXPECT_TRUE(IsFinite(2 * Gallon));
  EXPECT_FALSE(IsFinite((2 * Gallon) / l));
  EXPECT_FALSE(IsFinite((0 * Gallon) / l));
}

TEST_F(QuantitiesDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Quantity message;
    message.set_dimensions(0x7C00);
    message.set_magnitude(1.0);
    [[maybe_unused]] Speed const speed_of_light =
        Speed::ReadFromMessage(message);
  }, "representation.*dimensions");
}

TEST_F(QuantitiesTest, SerializationSuccess) {
  serialization::Quantity message;
  SpeedOfLight.WriteToMessage(&message);
  EXPECT_EQ(0x7C01, message.dimensions());
  EXPECT_EQ(299792458.0, message.magnitude());
  Speed const speed_of_light = Speed::ReadFromMessage(message);
  EXPECT_EQ(SpeedOfLight, speed_of_light);
}

// This check verifies that setting a log handler causes the protobuf library to
// report its errors using glog.  It doesn't have much too do with quantities,
// except that it's a convenient protobuf for this test.
TEST_F(QuantitiesDeathTest, SerializationLogHandler) {
  EXPECT_DEATH({
    google::protobuf::SetLogHandler(
        [](google::protobuf::LogLevel const level,
           char const* const filename,
           int const line,
           std::string const& message) {
          LOG_AT_LEVEL(level) << "[" << filename << ":" << line << "] "
                              << message;
        });
    serialization::Quantity message;
    message.set_magnitude(1.0);
    message.CheckInitialized();
  }, "missing required fields");
}

}  // namespace quantities
}  // namespace principia
