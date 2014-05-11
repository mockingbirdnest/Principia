#include "glog/logging.h"
#include "gtest/gtest.h"

#include "Geometry/Sign.hpp"

namespace principia {
namespace geometry {

class SignTest : public testing::Test {
 protected:
  Sign const positive_ = Sign(1);
  Sign const negative_ = Sign(-1);
};

TEST_F(SignTest, Integer) {
  LOG(INFO) << "JE SUIS LA!";
  Sign const positive(1);
  Sign const negative(-1);
  EXPECT_TRUE(positive.Positive());
  EXPECT_FALSE(positive.Negative());
  EXPECT_FALSE(negative.Positive());
  EXPECT_TRUE(negative.Negative());
}

TEST_F(SignTest, SignMultiplication) {
  Sign const positive(1);
  Sign const negative(-1);
  EXPECT_TRUE((positive * positive).Positive());
  EXPECT_TRUE((positive * negative).Negative());
  EXPECT_TRUE((negative * positive).Negative());
  EXPECT_TRUE((negative * negative).Positive());
}

TEST_F(SignTest, ScalarMultiplication) {
  Sign const positive(1);
  Sign const negative(-1);
  EXPECT_EQ(3, positive * 3);
  EXPECT_EQ(-3, positive * -3);
  EXPECT_EQ(-3, negative * 3);
  EXPECT_EQ(3, negative * -3);
}

}  // namespace geometry
}  // namespace principia
