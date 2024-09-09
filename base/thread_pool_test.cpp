#include "base/thread_pool.hpp"

#include <chrono>
#include <cstdint>
#include <future>
#include <thread>
#include <vector>

#include "absl/synchronization/mutex.h"
#include "absl/synchronization/notification.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

using namespace principia::base::_thread_pool;
using namespace std::chrono_literals;

class ThreadPoolTest : public ::testing::Test {
 protected:
  ThreadPoolTest() : pool_(std::thread::hardware_concurrency()) {
    LOG(ERROR) << "Concurrency is " << std::thread::hardware_concurrency();
  }

  ThreadPool<void> pool_;
};

// Check that execution occurs in parallel.  If things were sequential, the
// integers in `numbers` would be monotonically increasing.
TEST_F(ThreadPoolTest, ParallelExecution) {
#if defined(_DEBUG)
  constexpr int number_of_calls = 100'000;
#else
  constexpr int number_of_calls = 1'000'000;
#endif

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

TEST_F(ThreadPoolTest, TryAdd) {
  ThreadPool<void> pool(2);
  absl::Notification proceed1;
  absl::Notification proceed2;

  // Add two functions to the thread pool.  They get stuck and keep both threads
  // busy.
  pool.Add([&proceed1]() {  proceed1.WaitForNotification(); });
  pool.Add([&proceed2]() {  proceed2.WaitForNotification(); });

  // Try to add a call to the pool.  It doesn't work and the function doesn't
  // run.
  EXPECT_EQ(std::nullopt,
            pool.TryAdd([]() { LOG(FATAL) << "Should not run"; }));

  // Release one of the executing threads.  Now `TryAdd` should work and the
  // function should run.
  proceed1.Notify();
  std::this_thread::sleep_for(100ms);

  auto future = pool.TryAdd([]() { LOG(INFO) << "Hello, world!"; });
  future->wait();

  // Release the second thread.
  proceed2.Notify();
}

TEST_F(ThreadPoolTest, WaitUntilIdleFor) {
  ThreadPool<void> pool(2);
  absl::Notification proceed1;
  absl::Notification proceed2;

  // Add two functions to the thread pool.  They get stuck and keep both threads
  // busy.
  pool.Add([&proceed1]() {  proceed1.WaitForNotification(); });
  pool.Add([&proceed2]() {  proceed2.WaitForNotification(); });

  // Wait until a thread becomes idle.  That doesn't happen, so we timeout.
  EXPECT_FALSE(pool.WaitUntilIdleFor(absl::Milliseconds(100)));

  // Release one of the executing threads.  Now `WaitUntilIdleFor` should see
  // that there is an idle thread.
  proceed1.Notify();
  std::this_thread::sleep_for(100ms);

  EXPECT_TRUE(pool.WaitUntilIdleFor(absl::Hours(1)));

  // Release the second thread.
  proceed2.Notify();
}

}  // namespace base
}  // namespace principia
