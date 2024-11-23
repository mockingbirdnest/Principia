#include "base/push_pull_callback.hpp"

#include <random>
#include <utility>

#include "glog/logging.h"
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

TEST(PushPullCallback, Repeat) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> ordinate(-1.0, 1.0);

  double y1;
  double y2;

  auto f = [&y1, &y2](double const x) {
    // x = -1.0 -> y1; x = 1.0 -> y2.
    return (x + 1.0) * (y2 - y1) / 2.0 + y1;
  };

  auto task =
      [](std::function<double(double const x, void* const unused)> const& f)
      -> std::optional<double> {
    static constexpr std::int64_t n = 10000;
    double const f0 = f(0, nullptr);
    for (std::int64_t i = 1; i <= n; ++i) {
      double const fi = f(i / static_cast<double>(n), nullptr);
      if ((f0 > 0 && fi < 0) || (f0 < 0 && fi > 0)) {
        return i / static_cast<double>(n);
      }
    }
    return std::nullopt;
  };

  for (std::int64_t attempt = 0; attempt < 1000; ++ attempt) {
    y1 = ordinate(random);
    y2 = ordinate(random);

    auto* const executor =
        new PushPullExecutor<std::optional<double>, double, double, void*>(
            task);
    for (;;) {
      double x = -2;
      void* unused;
      bool const more = executor->callback().Pull(x, unused);
      if (!more) {
        auto const result = executor->get();
        if (result.has_value()) {
          LOG(ERROR) << attempt << ": " << result.value();
        } else {
          LOG(ERROR) << attempt << ": no result";
        }
        delete executor;
        break;
      }
      double y = f(x);
      executor->callback().Push(y);
    }
  }
}

}  // namespace base
}  // namespace principia
