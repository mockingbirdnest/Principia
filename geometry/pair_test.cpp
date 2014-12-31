#include "geometry/pair.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"

using principia::quantities::Action;
using principia::quantities::Winding;

namespace principia {
namespace geometry {

class PairTest : public testing::Test {
 protected:
  struct World;

  PairTest()
      : p1_(V1({4 * SIUnit<Action>(),
                5 * SIUnit<Action>(),
                6 * SIUnit<Action>()})),
        p2_(V2({16 * SIUnit<Winding>(),
                15 * SIUnit<Winding>(),
                14 * SIUnit<Winding>()})),
        v1_({1 * SIUnit<Action>(),
             2 * SIUnit<Action>(),
             3 * SIUnit<Action>()}),
        v2_({13 * SIUnit<Winding>(),
             12 * SIUnit<Winding>(),
             11 * SIUnit<Winding>()}),
        pp_(p1_, p2_),
        pv_(p1_, v2_),
        vp_(v1_, p2_),
        vv_(v1_, v2_) {}

  using V1 = Vector<Action, World>;
  using P1 = Point<V1>;
  using V2 = Vector<Winding, World>;
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

TEST_F(PairTest, Nothing) {
  //VP a(t1_, t2_);
  //P p2(t1_, t2_);
  //auto q = p1 - p2;
  //q = ((3.0 * (q - (-q))) * 2.0) / 5.0;
  //q *= 4.0;
  //p1 = p2 + q;
  //LOG(ERROR) << "p1 = " << p1 << "\nq = " << q;
}

TEST_F(PairTest, MemberAddition) {
  EXPECT_EQ(PP(P1(V1({5 * SIUnit<Action>(),
                      7 * SIUnit<Action>(),
                      9 * SIUnit<Action>()})),
               P2(V2({29 * SIUnit<Winding>(),
                      27 * SIUnit<Winding>(),
                      25 * SIUnit<Winding>()}))),
            pp_ + vv_);
  EXPECT_EQ(PV(P1(V1({5 * SIUnit<Action>(),
                      7 * SIUnit<Action>(),
                      9 * SIUnit<Action>()})),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            pv_ + vv_);
  EXPECT_EQ(VP(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               P2(V2({29 * SIUnit<Winding>(),
                      27 * SIUnit<Winding>(),
                      25 * SIUnit<Winding>()}))),
            vp_ + vv_);
  EXPECT_EQ(VV(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            vv_ + vv_);
}

}  // namespace geometry
}  // namespace principia
