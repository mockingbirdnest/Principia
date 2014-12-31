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

  using T1 = Vector<Action, World>;
  using T2 = Point<Vector<Winding, World>>;
  using P = Pair<T1, T2>;

  T1 t1_;
  T2 t2_;
};

TEST_F(PairTest, Nothing) {
  P p1(t1_, t2_);
  P p2(t1_, t2_);
  auto q = p1 - p2;
  p1 = p2 + q;
}

}  // namespace geometry
}  // namespace principia
