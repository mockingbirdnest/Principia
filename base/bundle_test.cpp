
#include "base/bundle.hpp"

#include <atomic>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "testing_utilities/matchers.hpp"

namespace principia {

using ::testing::Eq;

using namespace std::chrono_literals;

namespace base {

constexpr int workers = 8;

class BundleTest : public testing::Test {
 protected:
  Bundle bundle_;
};

using BundleDeathTest = BundleTest;

TEST_F(BundleDeathTest, AddAfterJoin) {
  EXPECT_DEATH({
      for (int i = 0; i < workers; ++i) {
        bundle_.Add([]() {
          std::this_thread::sleep_for(10ms);
          return Status::OK;
        });
      }
      bundle_.Join();
      bundle_.Add([]() {
        std::this_thread::sleep_for(10ms);
        return Status::OK;
      });
    },
    "!joining_");
}

TEST_F(BundleTest, MatrixVectorProduct) {
#if defined(_DEBUG)
  constexpr std::int64_t short_dimension = 100;
#else
  constexpr std::int64_t short_dimension = 1000;
#endif
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
        [&matrix, &vector, &product, i]() {
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

TEST_F(BundleTest, Deadline) {
  for (int i = 0; i < workers; ++i) {
    bundle_.Add([]() {
      std::this_thread::sleep_for(1000ms);
      return Status::OK;
    });
  }
  auto const status = bundle_.JoinWithin(10ms);
  EXPECT_THAT(status.error(), Eq(Error::DEADLINE_EXCEEDED));
  EXPECT_THAT(status.message(), Eq("bundle deadline exceeded"));
}

}  // namespace base
}  // namespace principia
