
#include "geometry/sign.hpp"

#include "glog/logging.h"
#include "gtest/gtest.h"

namespace principia {
namespace geometry {

class SignTest : public testing::Test {
 protected:
  Sign const positive_ = Sign::Positive();
  Sign const negative_ = Sign::Negative();
};

TEST_F(SignTest, Integer) {
  EXPECT_TRUE(positive_.is_positive());
  EXPECT_FALSE(positive_.is_negative());
  EXPECT_FALSE(negative_.is_positive());
  EXPECT_TRUE(negative_.is_negative());
}

TEST_F(SignTest, SignMultiplication) {
  EXPECT_TRUE((positive_ * positive_).is_positive());
  EXPECT_TRUE((positive_ * negative_).is_negative());
  EXPECT_TRUE((negative_ * positive_).is_negative());
  EXPECT_TRUE((negative_ * negative_).is_positive());
}

TEST_F(SignTest, ScalarMultiplication) {
  EXPECT_EQ(3, positive_ * 3);
  EXPECT_EQ(-3, positive_ * -3);
  EXPECT_EQ(-3, negative_ * 3);
  EXPECT_EQ(3, negative_ * -3);
}

TEST_F(SignTest, Serialization) {
  serialization::Sign message;
  Sign s = Sign::Positive();

  positive_.WriteToMessage(&message);
  EXPECT_FALSE(message.negative());
  s = Sign::ReadFromMessage(message);
  EXPECT_FALSE(s.is_negative());

  negative_.WriteToMessage(&message);
  EXPECT_TRUE(message.negative());
  s = Sign::ReadFromMessage(message);
  EXPECT_TRUE(s.is_negative());
}

}  // namespace geometry
}  // namespace principia
