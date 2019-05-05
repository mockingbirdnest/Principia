#include "physics/protector.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using geometry::Instant;
using quantities::si::Second;
using ::testing::MockFunction;

class ProtectorTest : public ::testing::Test {
 protected:
  MockFunction<void()> callback_;
  Protector protector_;
};

TEST_F(ProtectorTest, ImmediateExecution) {
  protector_.Protect(Instant() + 10 * Second);
  EXPECT_CALL(callback_, Call());
  CHECK(protector_.RunWhenUnprotected(Instant() + 5 * Second,
                                      callback_.AsStdFunction()));
}

TEST_F(ProtectorTest, DelayedExecution) {
  protector_.Protect(Instant() + 10 * Second);
  EXPECT_CALL(callback_, Call()).Times(0);
  CHECK(!protector_.RunWhenUnprotected(Instant() + 15 * Second,
                                       callback_.AsStdFunction()));

  protector_.Protect(Instant() + 20 * Second);

  EXPECT_CALL(callback_, Call()).Times(1);
  protector_.Unprotect(Instant() + 10 * Second);
}

}  // namespace physics
}  // namespace principia
