#include "ksp_plugin/pile_up.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin_test/mock_vessel.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using geometry::Displacement;
using geometry::Velocity;
using quantities::si::Kilo;
using quantities::si::Gram;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::Return;
using ::testing::ReturnRef;

class PileUpTest : public testing::Test {
 protected:
  MockVessel v1_;
  MockVessel v2_;
};

TEST_F(PileUpTest, ApparentDegreesOfFreedom) {
  Vector<Force, Barycentric> const null_intrinsic_force;

  Mass const mass1(1 * Kilo(Gram));
  Mass const mass2(2 * Kilo(Gram));

  DiscreteTrajectory<Barycentric> psychohistory1;
  DiscreteTrajectory<Barycentric> psychohistory2;
  psychohistory1.Append(
      astronomy::J2000,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({7 * Metre, 8 * Metre, 9 * Metre}),
          Velocity<Barycentric>({70 * Metre / Second,
                                 80 * Metre / Second,
                                 90 * Metre / Second})));
  psychohistory2.Append(
      astronomy::J2000,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({12 * Metre, 11 * Metre, 10 * Metre}),
          Velocity<Barycentric>({100 * Metre / Second,
                                 110 * Metre / Second,
                                 120 * Metre / Second})));

  EXPECT_CALL(v1_, psychohistory_is_authoritative())
      .WillRepeatedly(Return(true));
  EXPECT_CALL(v1_, intrinsic_force())
      .WillRepeatedly(ReturnRef(null_intrinsic_force));
  EXPECT_CALL(v1_, psychohistory()).WillRepeatedly(ReturnRef(psychohistory1));
  EXPECT_CALL(v1_, mass()).WillRepeatedly(ReturnRef(mass1));

  EXPECT_CALL(v2_, psychohistory_is_authoritative())
      .WillRepeatedly(Return(true));
  EXPECT_CALL(v2_, intrinsic_force())
      .WillRepeatedly(ReturnRef(null_intrinsic_force));
  EXPECT_CALL(v2_, psychohistory()).WillRepeatedly(ReturnRef(psychohistory2));
  EXPECT_CALL(v2_, mass()).WillRepeatedly(ReturnRef(mass2));

  PileUp pile_up({&v1_, &v2_});

  pile_up.SetVesselApparentDegreesOfFreedom(
      &v1_,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
          Velocity<Barycentric>({10 * Metre / Second,
                                 20 * Metre / Second,
                                 30 * Metre / Second})));
  pile_up.SetVesselApparentDegreesOfFreedom(
      &v2_,
      DegreesOfFreedom<Barycentric>(
          Barycentric::origin +
              Displacement<Barycentric>({6 * Metre, 5 * Metre, 4 * Metre}),
          Velocity<Barycentric>({60 * Metre / Second,
                                 50 * Metre / Second,
                                 40 * Metre / Second})));

  pile_up.UpdateVesselsInPileUp();

  EXPECT_THAT(
      pile_up.GetVesselActualDegreesOfFreedom(&v1_),
      Componentwise(
          AlmostEquals(
              Barycentric::origin +
                  Displacement<Barycentric>({7 * Metre, 8 * Metre, 9 * Metre}),
              1),
          AlmostEquals(Velocity<Barycentric>({(170.0 / 3.0) * Metre / Second,
                                              80 * Metre / Second,
                                              (310.0 / 3.0) * Metre / Second}),
                       1)));
  EXPECT_THAT(
      pile_up.GetVesselActualDegreesOfFreedom(&v2_),
      Componentwise(
          AlmostEquals(
              Barycentric::origin + Displacement<Barycentric>(
                                        {12 * Metre, 11 * Metre, 10 * Metre}),
              0),
          AlmostEquals(Velocity<Barycentric>({(320.0 / 3.0) * Metre / Second,
                                              110 * Metre / Second,
                                              (340.0 / 3.0) * Metre / Second}),
                       1)));
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
