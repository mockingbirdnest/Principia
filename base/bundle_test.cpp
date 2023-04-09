#include "base/bundle.hpp"

#include <atomic>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace base {

using ::testing::Eq;
using namespace principia::base::_bundle;
using namespace principia::testing_utilities::_matchers;
using namespace std::chrono_literals;

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
          return absl::OkStatus();
        });
      }
      bundle_.Join().IgnoreError();
      bundle_.Add([]() {
        std::this_thread::sleep_for(10ms);
        return absl::OkStatus();
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
        [i, &matrix, &product, &vector]() {
          product[i] = 0;
          for (int j = 0; j < long_dimension; ++j) {
            product[i] += matrix[i + short_dimension * j] * vector[j];
          }
          return absl::OkStatus();
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
      return absl::OkStatus();
    });
  }
  auto const status = bundle_.JoinWithin(10ms);
  EXPECT_THAT(status, StatusIs(absl::StatusCode::kDeadlineExceeded));
  EXPECT_THAT(status.message(), Eq("Bundle deadline exceeded"));
}

}  // namespace base
}  // namespace principia
