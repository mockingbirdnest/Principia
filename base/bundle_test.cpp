#include "base/bundle.hpp"

#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace principia {

using ::testing::Eq;

using namespace std::chrono_literals;  // NOLINT(build/namespaces)

namespace base {

class BundleTest : public testing::Test {
 protected:
  BundleTest() : bundle_(/*workers=*/8) {}

  Bundle bundle_;
};

TEST_F(BundleTest, MatrixVectorProduct) {
  constexpr std::int64_t short_dimension = 100;
  constexpr std::int64_t long_dimension = 100000;
  std::vector<std::int64_t> matrix(short_dimension * long_dimension, 1);
  std::vector<std::int64_t> vector(long_dimension);
  std::vector<std::int64_t> product(short_dimension);
  for (int i = 0; i < short_dimension; ++i) {
    matrix[i + short_dimension * i] = 2;
  }
  for (int j = 0; j < long_dimension; ++j) {
    vector[j] = j;
  }
  for (int i = 0; i < short_dimension; ++i) {
    bundle_.Add(
        [&matrix, &vector, &product, i, short_dimension, long_dimension]() {
          product[i] = 0;
          for (int j = 0; j < long_dimension; ++j) {
            product[i] += matrix[i + short_dimension * j] * vector[j];
          }
          return Status::OK;
        });
  }
  EXPECT_OK(bundle_.Join());
  for (int i = 0; i < short_dimension; ++i) {
    EXPECT_THAT(product[i],
                Eq(long_dimension * (long_dimension - 1) / 2 + i));
  }
}

TEST_F(BundleTest, Abort) {
  int i = 0;
  auto const wait = [this]() {
    while (!AbortRequested()) {
      std::this_thread::sleep_for(10ms);
    }
    return Status(Error::ABORTED, "");
  };
  auto const wait_with_message = [this]() {
    while (!AbortRequested()) {
      std::this_thread::sleep_for(10ms);
    }
    return Status(Error::ABORTED, "aborted on request");
  };
  auto const increment_i = [&i]() {
    ++i;
    return Status::OK;
  };
  std::vector<Bundle::TaskHandle> waiters;
  Bundle::TaskHandle waiter;
  for (int i = 0; i < 7; ++i) {
    waiters.emplace_back(bundle_.Add(wait));
  }
  waiter = bundle_.Add(wait_with_message);
  // All threads busy.
  Bundle::TaskHandle not_scheduled = bundle_.Add(increment_i);
  std::this_thread::sleep_for(10ms);
  bundle_.Abort(not_scheduled);
  EXPECT_THAT(i, Eq(0));
  bundle_.Abort(waiter);
  // Now a thread is free.
  Bundle::TaskHandle incrementer = bundle_.Add(increment_i);
  EXPECT_OK(bundle_.JoinTask(std::move(incrementer)));
  EXPECT_THAT(i, Eq(1));
  // But only one.
  waiter = bundle_.Add(wait);
  not_scheduled = bundle_.Add(increment_i);
  std::this_thread::sleep_for(10ms);
  bundle_.Abort(not_scheduled);
  EXPECT_THAT(i, Eq(1));
  bundle_.Abort(waiter);
  for (auto const& w : waiters) {
    bundle_.Abort(w);
  }
  auto status = bundle_.Join();
  EXPECT_TRUE(status.error() == Error::ABORTED);
  EXPECT_THAT(status.message(), Eq("aborted on request"));
}

}  // namespace base
}  // namespace principia
