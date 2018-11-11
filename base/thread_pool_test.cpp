
#include "base/thread_pool.hpp"

#include <vector>

#include "absl/synchronization/mutex.h"
#include "glog/logging.h"
#include "gmock/gmock.h"

namespace principia {
namespace base {

class ThreadPoolTest : public ::testing::Test {
 protected:
  ThreadPoolTest() : pool_(std::thread::hardware_concurrency()) {
    LOG(ERROR) << "Concurrency is " << std::thread::hardware_concurrency();
  }

  ThreadPool<void> pool_;
};

// Check that execution occurs in parallel.  If things were sequential, the
// integers in |numbers| would be monotonically increasing.
TEST_F(ThreadPoolTest, ParallelExecution) {
  static constexpr int number_of_calls = 1'000'000;

  absl::Mutex lock;
  std::vector<std::int64_t> numbers;
  std::vector<std::future<void>> futures;
  for (std::int64_t i = 0; i < number_of_calls; ++i) {
    futures.push_back(pool_.Add([i, &lock, &numbers]() {
      absl::MutexLock l(&lock);
      numbers.push_back(i);
    }));
  }

  for (auto const& future : futures) {
    future.wait();
  }

  bool monotonically_increasing = true;
  for (std::int64_t i = 1; i < numbers.size(); ++i) {
    if (numbers[i] < numbers[i - 1]) {
      monotonically_increasing = false;
    }
  }
  EXPECT_FALSE(monotonically_increasing);
}

}  // namespace base
}  // namespace principia
