#include "base/push_pull_callback.hpp"

#include <utility>

#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace principia::base::_push_pull_callback;

TEST(PushPullCallback, Test) {
  // A task that does some computation, including running its callback `f` with
  // various arguments.
  auto task = [](std::function<int(int left, int right)> const& f) {
    double const a = f(2, 4);
    double const b = f(3, 5);
    return a - b;
  };

  PushPullExecutor<double, int, int, int> executor(std::move(task));
  auto& callback = executor.callback();

  int left;
  int right;

  EXPECT_TRUE(callback.Pull(left, right));
  EXPECT_EQ(2, left);
  EXPECT_EQ(4, right);
  callback.Push(left + right);

  EXPECT_TRUE(callback.Pull(left, right));
  EXPECT_EQ(3, left);
  EXPECT_EQ(5, right);
  callback.Push(left - right);

  EXPECT_FALSE(callback.Pull(left, right));
  EXPECT_EQ(8, executor.get());
}

TEST(PushPullCallback, 4136) {
  auto task =
      [](std::function<double(int const x)> const& f)
      -> std::optional<double> {
    return std::nullopt;
  };

  for (std::int64_t attempt = 0; attempt < 100'000; ++attempt) {
    auto* const executor =
        new PushPullExecutor<std::optional<double>, double, int>(task);
    for (;;) {
      int x;
      bool const more = executor->callback().Pull(x);
      if (!more) {
        auto const result = executor->get();
        delete executor;
        break;
      }
      executor->callback().Push(1.0);
    }
  }
}

}  // namespace base
}  // namespace principia
