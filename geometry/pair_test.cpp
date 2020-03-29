
#include "geometry/pair.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

// Define this to check that the illegal operations are actually rejected.
#undef CHECK_ILLEGAL

namespace principia {

using quantities::Action;
using quantities::Amount;
using quantities::Angle;
using quantities::CatalyticActivity;
using quantities::Energy;
using quantities::Entropy;
using quantities::Product;
using quantities::Quantity;
using quantities::SolidAngle;
using quantities::Time;
namespace si = quantities::si;

namespace geometry {

class PairTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Universe = Frame<enum class UniverseTag>;

  PairTest()
      : p1_(P1() + V1({4 * si::Unit<Action>,
                       5 * si::Unit<Action>,
                       6 * si::Unit<Action>})),
        p2_(P2() + V2({16 * si::Unit<Amount>,
                       15 * si::Unit<Amount>,
                       14 * si::Unit<Amount>})),
        v1_({1 * si::Unit<Action>,
             2 * si::Unit<Action>,
             3 * si::Unit<Action>}),
        v2_({13 * si::Unit<Amount>,
             12 * si::Unit<Amount>,
             11 * si::Unit<Amount>}),
        pp_(p1_, p2_),
        pv_(p1_, v2_),
        vp_(v1_, p2_),
        vv_(v1_, v2_) {}

  using V1 = Vector<Action, World>;
  using P1 = Point<V1>;
  using V2 = Vector<Amount, World>;
  using P2 = Point<V2>;
  using PP = Pair<P1, P2>;
  using PV = Pair<P1, V2>;
  using VP = Pair<V1, P2>;
  using VV = Pair<V1, V2>;

  P1 p1_;
  P2 p2_;
  V1 v1_;
  V2 v2_;
  PP pp_;
  PV pv_;
  VP vp_;
  VV vv_;
};

using PairDeathTest = PairTest;

TEST_F(PairTest, MemberAddition) {
  EXPECT_EQ(PP(P1() + V1({5 * si::Unit<Action>,
                          7 * si::Unit<Action>,
                          9 * si::Unit<Action>}),
               P2() + V2({29 * si::Unit<Amount>,
                          27 * si::Unit<Amount>,
                          25 * si::Unit<Amount>})),
            pp_ + vv_);
  EXPECT_EQ(PV(P1() + V1({5 * si::Unit<Action>,
                          7 * si::Unit<Action>,
                          9 * si::Unit<Action>}),
               V2({26 * si::Unit<Amount>,
                   24 * si::Unit<Amount>,
                   22 * si::Unit<Amount>})),
            pv_ + vv_);
  EXPECT_EQ(VP(V1({2 * si::Unit<Action>,
                   4 * si::Unit<Action>,
                   6 * si::Unit<Action>}),
               P2() + V2({29 * si::Unit<Amount>,
                          27 * si::Unit<Amount>,
                          25 * si::Unit<Amount>})),
            vp_ + vv_);
  EXPECT_EQ(VV(V1({2 * si::Unit<Action>,
                   4 * si::Unit<Action>,
                   6 * si::Unit<Action>}),
               V2({26 * si::Unit<Amount>,
                   24 * si::Unit<Amount>,
                   22 * si::Unit<Amount>})),
            vv_ + vv_);
}

TEST_F(PairTest, MemberSubtraction) {
  EXPECT_EQ(PP(P1() + V1({3 * si::Unit<Action>,
                          3 * si::Unit<Action>,
                          3 * si::Unit<Action>}),
               P2() + V2({3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>})),
            pp_ - vv_);
  EXPECT_EQ(PV(P1() + V1({3 * si::Unit<Action>,
                          2 * si::Unit<Action>,
                          3 * si::Unit<Action>}) +
                      // A convoluted way of writing {3, 3, 3}, to circumvent a
                      // frightening bug in VS2017.
                      V1({0 * si::Unit<Action>,
                          1 * si::Unit<Action>,
                          0 * si::Unit<Action> }),
               V2({0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>})),
            pv_ - vv_);
  EXPECT_EQ(VP(V1({0 * si::Unit<Action>,
                   0 * si::Unit<Action>,
                   0 * si::Unit<Action>}),
               P2() + V2({3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>})),
            vp_ - vv_);
  EXPECT_EQ(VV(V1({0 * si::Unit<Action>,
                   0 * si::Unit<Action>,
                   0 * si::Unit<Action>}),
               V2({0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>})),
            vv_ - vv_);
}

TEST_F(PairTest, MemberAdditionTo) {
  pp_ += vv_;
  EXPECT_EQ(PP(P1() + V1({5 * si::Unit<Action>,
                          7 * si::Unit<Action>,
                          9 * si::Unit<Action>}),
               P2() + V2({29 * si::Unit<Amount>,
                          27 * si::Unit<Amount>,
                          25 * si::Unit<Amount>})),
            pp_);
  pv_ += vv_;
  EXPECT_EQ(PV(P1() + V1({5 * si::Unit<Action>,
                          7 * si::Unit<Action>,
                          9 * si::Unit<Action>}),
               V2({26 * si::Unit<Amount>,
                   24 * si::Unit<Amount>,
                   22 * si::Unit<Amount>})),
            pv_);
  vp_ += vv_;
  EXPECT_EQ(VP(V1({2 * si::Unit<Action>,
                   4 * si::Unit<Action>,
                   6 * si::Unit<Action>}),
               P2() + V2({29 * si::Unit<Amount>,
                          27 * si::Unit<Amount>,
                          25 * si::Unit<Amount>})),
            vp_);
  vv_ += vv_;
  EXPECT_EQ(VV(V1({2 * si::Unit<Action>,
                   4 * si::Unit<Action>,
                   6 * si::Unit<Action>}),
               V2({26 * si::Unit<Amount>,
                   24 * si::Unit<Amount>,
                   22 * si::Unit<Amount>})),
            vv_);
}

TEST_F(PairTest, MemberSubtractionFrom) {
  pp_ -= vv_;
  EXPECT_EQ(PP(P1() + V1({3 * si::Unit<Action>,
                          3 * si::Unit<Action>,
                          3 * si::Unit<Action>}),
               P2() + V2({3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>})),
            pp_);
  pv_ -= vv_;
  EXPECT_EQ(PV(P1() + V1({3 * si::Unit<Action>,
                          3 * si::Unit<Action>,
                          3 * si::Unit<Action>}),
               V2({0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>})),
            pv_);
  vp_ -= vv_;
  EXPECT_EQ(VP(V1({0 * si::Unit<Action>,
                   0 * si::Unit<Action>,
                   0 * si::Unit<Action>}),
               P2() + V2({3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>,
                          3 * si::Unit<Amount>})),
            vp_);
  vv_ -= vv_;
  EXPECT_EQ(VV(V1({0 * si::Unit<Action>,
                   0 * si::Unit<Action>,
                   0 * si::Unit<Action>}),
               V2({0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>,
                   0 * si::Unit<Amount>})),
            vv_);
}

TEST_F(PairTest, MemberEquality) {
  EXPECT_TRUE(pp_ == PP(P1() + V1({4 * si::Unit<Action>,
                                   5 * si::Unit<Action>,
                                   6 * si::Unit<Action>}),
                        P2() + V2({16 * si::Unit<Amount>,
                                   15 * si::Unit<Amount>,
                                   14 * si::Unit<Amount>})));
  EXPECT_FALSE(pp_ == PP(P1() + V1({4 * si::Unit<Action>,
                                    5 * si::Unit<Action>,
                                    6 * si::Unit<Action>}),
                         P2() + V2({26 * si::Unit<Amount>,
                                15 * si::Unit<Amount>,
                                14 * si::Unit<Amount>})));
  EXPECT_TRUE(pv_ == PV(P1() + V1({4 * si::Unit<Action>,
                                   5 * si::Unit<Action>,
                                   6 * si::Unit<Action>}),
                        V2({13 * si::Unit<Amount>,
                            12 * si::Unit<Amount>,
                            11 * si::Unit<Amount>})));
  EXPECT_FALSE(pv_ == PV(P1() + V1({4 * si::Unit<Action>,
                                    15 * si::Unit<Action>,
                                    6 * si::Unit<Action>}),
                         V2({13 * si::Unit<Amount>,
                             12 * si::Unit<Amount>,
                             11 * si::Unit<Amount>})));
  EXPECT_TRUE(vp_ == VP(V1({1 * si::Unit<Action>,
                            2 * si::Unit<Action>,
                            3 * si::Unit<Action>}),
                        P2() + V2({16 * si::Unit<Amount>,
                               15 * si::Unit<Amount>,
                               14 * si::Unit<Amount>})));
  EXPECT_FALSE(vp_ == VP(V1({1 * si::Unit<Action>,
                             2 * si::Unit<Action>,
                             13 * si::Unit<Action>}),
                         P2() + V2({16 * si::Unit<Amount>,
                                15 * si::Unit<Amount>,
                                14 * si::Unit<Amount>})));
  EXPECT_TRUE(vv_ == VV(V1({1 * si::Unit<Action>,
                            2 * si::Unit<Action>,
                            3 * si::Unit<Action>}),
                        V2({13 * si::Unit<Amount>,
                            12 * si::Unit<Amount>,
                            11 * si::Unit<Amount>})));
  EXPECT_FALSE(vv_ == VV(V1({1 * si::Unit<Action>,
                             2 * si::Unit<Action>,
                             3 * si::Unit<Action>}),
                         V2({13 * si::Unit<Amount>,
                             12 * si::Unit<Amount>,
                             1 * si::Unit<Amount>})));
}

TEST_F(PairTest, MemberInequality) {
  EXPECT_FALSE(pp_ != PP(P1() + V1({4 * si::Unit<Action>,
                                    5 * si::Unit<Action>,
                                    6 * si::Unit<Action>}),
                         P2() + V2({16 * si::Unit<Amount>,
                                    15 * si::Unit<Amount>,
                                    14 * si::Unit<Amount>})));
  EXPECT_TRUE(pp_ != PP(P1() + V1({4 * si::Unit<Action>,
                                   5 * si::Unit<Action>,
                                   6 * si::Unit<Action>}),
                        P2() + V2({26 * si::Unit<Amount>,
                                   15 * si::Unit<Amount>,
                                   14 * si::Unit<Amount>})));
  EXPECT_FALSE(pv_ != PV(P1() + V1({4 * si::Unit<Action>,
                                    5 * si::Unit<Action>,
                                    6 * si::Unit<Action>}),
                         V2({13 * si::Unit<Amount>,
                             12 * si::Unit<Amount>,
                             11 * si::Unit<Amount>})));
  EXPECT_TRUE(pv_ != PV(P1() + V1({4 * si::Unit<Action>,
                                   15 * si::Unit<Action>,
                                   6 * si::Unit<Action>}),
                        V2({13 * si::Unit<Amount>,
                            12 * si::Unit<Amount>,
                            11 * si::Unit<Amount>})));
  EXPECT_FALSE(vp_ != VP(V1({1 * si::Unit<Action>,
                             2 * si::Unit<Action>,
                             3 * si::Unit<Action>}),
                         P2() + V2({16 * si::Unit<Amount>,
                                15 * si::Unit<Amount>,
                                14 * si::Unit<Amount>})));
  EXPECT_TRUE(vp_ != VP(V1({1 * si::Unit<Action>,
                            2 * si::Unit<Action>,
                            13 * si::Unit<Action>}),
                        P2() + V2({16 * si::Unit<Amount>,
                               15 * si::Unit<Amount>,
                               14 * si::Unit<Amount>})));
  EXPECT_FALSE(vv_ != VV(V1({1 * si::Unit<Action>,
                             2 * si::Unit<Action>,
                             3 * si::Unit<Action>}),
                         V2({13 * si::Unit<Amount>,
                             12 * si::Unit<Amount>,
                             11 * si::Unit<Amount>})));
  EXPECT_TRUE(vv_ != VV(V1({1 * si::Unit<Action>,
                            2 * si::Unit<Action>,
                            3 * si::Unit<Action>}),
                        V2({13 * si::Unit<Amount>,
                            12 * si::Unit<Amount>,
                            1 * si::Unit<Amount>})));
}

TEST_F(PairTest, AffineSubtraction) {
  PP const pp = pp_ + vv_;
  EXPECT_EQ(vv_, pp - pp_);
  PV const pv = pv_ + vv_;
  EXPECT_EQ(vv_, pv - pv_);
  VP const vp = vp_ + vv_;
  EXPECT_EQ(vv_, vp - vp_);
  // No test for VV, that would be a vector subtraction.
}

TEST_F(PairTest, UnaryPlus) {
  VV const vv = +vv_;
  EXPECT_EQ(vv_, vv);
#ifdef CHECK_ILLEGAL
  auto const pp = +pp_;
  auto const pv = +pv_;
  auto const vp = +vp_;
#endif
}

TEST_F(PairTest, UnaryMinus) {
  VV const vv = -vv_;
  EXPECT_EQ(VV(V1({-1 * si::Unit<Action>,
                   -2 * si::Unit<Action>,
                   -3 * si::Unit<Action>}),
               V2({-13 * si::Unit<Amount>,
                   -12 * si::Unit<Amount>,
                   -11 * si::Unit<Amount>})),
            vv);
#ifdef CHECK_ILLEGAL
  auto const pp = -pp_;
  auto const pv = -pv_;
  auto const vp = -vp_;
#endif
}

TEST_F(PairTest, LeftMultiplication) {
  EXPECT_EQ(VV(V1({3 * si::Unit<Action>,
                   6 * si::Unit<Action>,
                   9 * si::Unit<Action>}),
               V2({39 * si::Unit<Amount>,
                   36 * si::Unit<Amount>,
                   33 * si::Unit<Amount>})),
            3 * vv_);
#ifdef CHECK_ILLEGAL
  auto const pp = 3 * pp_;
  auto const pv = 3 * pv_;
  auto const vp = 3 * vp_;
#endif
}

TEST_F(PairTest, RightMultiplication) {
  using Whatever1 = Product<Action, SolidAngle>;
  using Whatever2 = Product<Amount, SolidAngle>;
  using VWhatever1 = Vector<Whatever1, World>;
  using VWhatever2 = Vector<Whatever2, World>;
  using Pear = Pair<VWhatever1, VWhatever2>;
  EXPECT_EQ(Pear(VWhatever1({1.5 * si::Unit<Whatever1>,
                             3.0 * si::Unit<Whatever1>,
                             4.5 * si::Unit<Whatever1>}),
                 VWhatever2({19.5 * si::Unit<Whatever2>,
                             18.0 * si::Unit<Whatever2>,
                             16.5 * si::Unit<Whatever2>})),
            vv_ * (1.5 * si::Unit<SolidAngle>));
#ifdef CHECK_ILLEGAL
  auto const pp = pp_ * (1.5 * si::Unit<SolidAngle>);
  auto const pv = pv_ * (1.5 * si::Unit<SolidAngle>);
  auto const vp = vp_ * (1.5 * si::Unit<SolidAngle>);
#endif
}

TEST_F(PairTest, RightDivision) {
  using VEnergy = Vector<Energy, World>;
  using VCatalyticActivity = Vector<CatalyticActivity, World>;
  using Pear = Pair<VEnergy, VCatalyticActivity>;
  EXPECT_EQ(Pear(VEnergy({0.5 * si::Unit<Energy>,
                          1.0 * si::Unit<Energy>,
                          1.5 * si::Unit<Energy>}),
                 VCatalyticActivity({6.5 * si::Unit<CatalyticActivity>,
                                     6.0 * si::Unit<CatalyticActivity>,
                                     5.5 * si::Unit<CatalyticActivity>})),
            vv_ / (2 * si::Unit<Time>));
#ifdef CHECK_ILLEGAL
  auto const pp = pp_ / (2 * si::Unit<Time>);
  auto const pv = pv_ / (2 * si::Unit<Time>);
  auto const vp = vp_ / (2 * si::Unit<Time>);
#endif
}

TEST_F(PairTest, MultiplicationBy) {
  vv_ *= 4;
  EXPECT_EQ(VV(V1({4 * si::Unit<Action>,
                   8 * si::Unit<Action>,
                   12 * si::Unit<Action>}),
               V2({52 * si::Unit<Amount>,
                   48 * si::Unit<Amount>,
                   44 * si::Unit<Amount>})),
            vv_);
#ifdef CHECK_ILLEGAL
  pp_ *= 4;
  pv_ *= 4;
  vp_ *= 4;
#endif
}

TEST_F(PairTest, DivisionBy) {
  vv_ /= 0.25;
  EXPECT_EQ(VV(V1({4 * si::Unit<Action>,
                   8 * si::Unit<Action>,
                   12 * si::Unit<Action>}),
               V2({52 * si::Unit<Amount>,
                   48 * si::Unit<Amount>,
                   44 * si::Unit<Amount>})),
            vv_);
#ifdef CHECK_ILLEGAL
  pp_ /= 0.25;
  pv_ /= 0.25;
  vp_ /= 0.25;
#endif
}

TEST_F(PairTest, Streaming) {
  LOG(ERROR) << "pp_ = " << pp_;
  LOG(ERROR) << "pv_ = " << pv_;
  LOG(ERROR) << "vp_ = " << vp_;
  LOG(ERROR) << "vv_ = " << vv_;
}

TEST_F(PairDeathTest, SerializationError) {
  serialization::Pair message;
  pv_.WriteToMessage(&message);
  EXPECT_DEATH({
    PP const pp = PP::ReadFromMessage(message);
  }, "has_point");
  EXPECT_DEATH({
    VP const vp = VP::ReadFromMessage(message);
  }, "has_multivector");
  EXPECT_DEATH({
    VV const vv = VV::ReadFromMessage(message);
  }, "has_multivector");
}

TEST_F(PairTest, SerializationSuccess) {
  serialization::Pair message;

  pp_.WriteToMessage(&message);
  EXPECT_TRUE(message.t1().has_point());
  EXPECT_FALSE(message.t1().has_multivector());
  EXPECT_TRUE(message.t2().has_point());
  EXPECT_FALSE(message.t2().has_multivector());
  PP const pp = PP::ReadFromMessage(message);
  EXPECT_EQ(pp_, pp);

  pv_.WriteToMessage(&message);
  EXPECT_TRUE(message.t1().has_point());
  EXPECT_FALSE(message.t1().has_multivector());
  EXPECT_FALSE(message.t2().has_point());
  EXPECT_TRUE(message.t2().has_multivector());
  PV const pv = PV::ReadFromMessage(message);
  EXPECT_EQ(pv_, pv);

  vp_.WriteToMessage(&message);
  EXPECT_FALSE(message.t1().has_point());
  EXPECT_TRUE(message.t1().has_multivector());
  EXPECT_TRUE(message.t2().has_point());
  EXPECT_FALSE(message.t2().has_multivector());
  VP const vp = VP::ReadFromMessage(message);
  EXPECT_EQ(vp_, vp);

  vv_.WriteToMessage(&message);
  EXPECT_FALSE(message.t1().has_point());
  EXPECT_TRUE(message.t1().has_multivector());
  EXPECT_FALSE(message.t2().has_point());
  EXPECT_TRUE(message.t2().has_multivector());
  VV const vv = VV::ReadFromMessage(message);
  EXPECT_EQ(vv_, vv);
}

TEST_F(PairDeathTest, BarycentreCalculatorError) {
  using PPBarycentreCalculator = BarycentreCalculator<PP, Entropy>;
  using PVBarycentreCalculator = BarycentreCalculator<PV, Entropy>;
  using VPBarycentreCalculator = BarycentreCalculator<VP, Entropy>;
  using VVBarycentreCalculator = BarycentreCalculator<VV, Entropy>;
  EXPECT_DEATH({
    PPBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    PVBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    VPBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    VVBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

// Because we can't access individual members this test doesn't fully exercise
// the computations, so we'll redo some testing for DegreesOfFreedom.
TEST_F(PairTest, BarycentreCalculatorSuccess) {
  {
    BarycentreCalculator<PP, double> calculator;
    calculator.Add(pp_, 3);
    PP barycentre = calculator.Get();
    EXPECT_EQ(pp_, barycentre);
    calculator.Add(pp_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(pp_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(pp_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(pp_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<PV, double> calculator;
    calculator.Add(pv_, 3);
    PV barycentre = calculator.Get();
    EXPECT_EQ(pv_, barycentre);
    calculator.Add(pv_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(pv_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(pv_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(pv_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<VP, double> calculator;
    calculator.Add(vp_, 3);
    VP barycentre = calculator.Get();
    EXPECT_EQ(vp_, barycentre);
    calculator.Add(vp_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(vp_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(vp_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(vp_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<VV, double> calculator;
    calculator.Add(vv_, 3);
    VV barycentre = calculator.Get();
    EXPECT_EQ(vv_, barycentre);
    calculator.Add(vv_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(13.0 / 8.0 * vv_, barycentre);
    calculator.Add(-2 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(-3.0 / 16.0 * vv_, barycentre);
  }
}

TEST_F(PairTest, Mappable) {
  using VV1U = Vector<Action, Universe>;
  using VV2U = Vector<Amount, Universe>;
  using VVU = Pair<VV1U, VV2U>;
  VVU const vv = Identity<World, Universe>()(vv_);
  EXPECT_EQ(VVU(VV1U({1 * si::Unit<Action>,
                      2 * si::Unit<Action>,
                      3 * si::Unit<Action>}),
                VV2U({13 * si::Unit<Amount>,
                      12 * si::Unit<Amount>,
                      11 * si::Unit<Amount>})),
            vv);

#ifdef CHECK_ILLEGAL
  auto const pp = Identity<World, Universe>()(pp_);
  auto const pv = Identity<World, Universe>()(pv_);
  auto const vp = Identity<World, Universe>()(vp_);
#endif
}

}  // namespace geometry
}  // namespace principia
