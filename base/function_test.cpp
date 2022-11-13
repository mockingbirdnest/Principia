#include "base/function.hpp"

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::Eq;

namespace principia {
namespace base {

class FunctionTest : public testing::Test {};

TEST_F(FunctionTest, MovableFunction) {
  int observe_i = 0;
  function<void()> λ = [i = std::make_unique<int>(), &observe_i](){
    ++*i;
    observe_i = *i;
  };
  function<void()> f = std::move(λ);
  f();
  EXPECT_THAT(observe_i, Eq(1));
  f();
  EXPECT_THAT(observe_i, Eq(2));

  function<std::unique_ptr<int>(std::unique_ptr<int>)> μ =
      [i = std::make_unique<int>()](std::unique_ptr<int> p) {
    ++*i;
    *p += *i;
    return p;
  };
  EXPECT_THAT(*μ(std::make_unique<int>(3)), Eq(4));
  EXPECT_THAT(*μ(std::make_unique<int>(0)), Eq(2));
}

}  // namespace base
}  // namespace principia
