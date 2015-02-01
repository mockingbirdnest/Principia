#include "geometry/sign.hpp"

#include "glog/logging.h"
#include "gtest/gtest.h"

namespace principia {
namespace geometry {

class SignTest : public testing::Test {
 protected:
  Sign const positive_ = Sign(1);
  Sign const negative_ = Sign(-1);
};

TEST_F(SignTest, Integer) {
  EXPECT_TRUE(positive_.Positive());
  EXPECT_FALSE(positive_.Negative());
  EXPECT_FALSE(negative_.Positive());
  EXPECT_TRUE(negative_.Negative());
}

TEST_F(SignTest, SignMultiplication) {
  EXPECT_TRUE((positive_ * positive_).Positive());
  EXPECT_TRUE((positive_ * negative_).Negative());
  EXPECT_TRUE((negative_ * positive_).Negative());
  EXPECT_TRUE((negative_ * negative_).Positive());
}

TEST_F(SignTest, ScalarMultiplication) {
  EXPECT_EQ(3, positive_ * 3);
  EXPECT_EQ(-3, positive_ * -3);
  EXPECT_EQ(-3, negative_ * 3);
  EXPECT_EQ(3, negative_ * -3);
}

TEST_F(SignTest, Serialization) {
  serialization::Sign message;
  Sign s(1);

  positive_.WriteToMessage(&message);
  EXPECT_FALSE(message.negative());
  s = Sign::ReadFromMessage(message);
  EXPECT_FALSE(s.Negative());

  negative_.WriteToMessage(&message);
  EXPECT_TRUE(message.negative());
  s = Sign::ReadFromMessage(message);
  EXPECT_TRUE(s.Negative());
}

}  // namespace geometry
}  // namespace principia
