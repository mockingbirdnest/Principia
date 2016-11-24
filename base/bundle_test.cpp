#include "base/bundle.hpp"

#include <atomic>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace principia {

using ::testing::Eq;

using namespace std::chrono_literals;  // NOLINT(build/namespaces)

namespace base {

constexpr int workers = 8;

class BundleTest : public testing::Test {
 protected:
  BundleTest() : bundle_(workers) {}

  Status Wait() {
    ++waiters_activated_;
    while (!AbortRequested()) {
      std::this_thread::sleep_for(10ms);
    }
    ++waiters_terminated_;
    return Status(Error::ABORTED, "");
  };

  Bundle bundle_;
  std::atomic<int> waiters_activated_ = 0;
  std::atomic<int> waiters_terminated_ = 0;
  Bundle::Task const wait_ = std::bind(&BundleTest::Wait, this);
  Bundle::Task const cancel_ = [] { return Status::CANCELLED; };
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
  for (int i = 0; i < workers - 1; ++i) {
    bundle_.Add(wait_);
  }
  while (waiters_activated_ < workers - 1) {
    std::this_thread::sleep_for(10ms);
  }
  EXPECT_THAT(waiters_terminated_, Eq(0));
  bundle_.Add(cancel_);
  // This one will never start.
  bundle_.Add(wait_);
  EXPECT_THAT(bundle_.Join(), Eq(Status::CANCELLED));
  EXPECT_THAT(waiters_activated_, Eq(workers - 1));
  EXPECT_THAT(waiters_terminated_, Eq(workers - 1));
}

TEST_F(BundleTest, NestedAbort) {
  constexpr int workers_per_dependent_bundle = 2;
  auto const own_a_bundle = [this, workers_per_dependent_bundle] {
    Bundle dependent(workers_per_dependent_bundle);
    // |workers_per_dependent_bundle| tasks will run, and an additional one will
    // be queued.
    for (int i = 0; i < workers_per_dependent_bundle + 1; ++i) {
      dependent.Add(wait_);
    }
    auto status = dependent.Join();
    EXPECT_OK(status);
    return status;
  };

  for (int i = 0; i < workers - 1; ++i) {
    bundle_.Add(own_a_bundle);
  }
  while (waiters_activated_ < (workers - 1) * workers_per_dependent_bundle) {
    std::this_thread::sleep_for(10ms);
  }
  EXPECT_THAT(waiters_terminated_, Eq(0));
  bundle_.Add(cancel_);
  EXPECT_THAT(bundle_.Join(), Eq(Status::CANCELLED));
  EXPECT_THAT(waiters_activated_,
              Eq((workers - 1) * workers_per_dependent_bundle));
  EXPECT_THAT(waiters_terminated_,
              Eq((workers - 1) * workers_per_dependent_bundle));
}

}  // namespace base
}  // namespace principia
