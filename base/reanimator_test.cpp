#include "base/reanimator.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <thread>
#include <vector>

#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"
#include "absl/synchronization/notification.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/matchers.hpp"

namespace principia {
namespace base {

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;
using namespace principia::base::_reanimator;
using namespace principia::testing_utilities::_matchers;
using namespace std::chrono_literals;

class ReanimatorTest : public ::testing::Test {
 protected:
  using ToyReanimator = Reanimator<int>;
};

TEST_F(ReanimatorTest, RunGuaranteed) {
  // No need for locking, at most one action is running at any point in time.
  std::vector<int> processed;
  ToyReanimator reanimator([&processed](int const key) {
    processed.push_back(key);
    return absl::OkStatus();
  });

  auto const handle1 = reanimator.RunGuaranteed(1);
  auto const handle2 = reanimator.RunGuaranteed(2);
  auto const handle3 = reanimator.RunGuaranteed(3);

  // Starting the reanimator *after* queueing the actions ensures that they are
  // executed in order.
  reanimator.Start();

  // Wait for the actions in a order different from their actual execution
  // order.  When the action 1 completes, surely all actions have completed.
  EXPECT_OK(reanimator.Wait(handle2));
  EXPECT_OK(reanimator.Wait(handle1));

  // It's fine to wait on a stopped reanimator.
  reanimator.Stop();
  EXPECT_OK(reanimator.Wait(handle3));

  // The runs execute in decreasing order of key.
  EXPECT_THAT(processed, ElementsAre(3, 2, 1));
}

TEST_F(ReanimatorTest, RunBestEffort) {
  std::vector<int> processed;
  ToyReanimator reanimator([&processed](int const& key) {
    processed.push_back(key);
    return absl::OkStatus();
  });

  reanimator.Start();

  reanimator.RunBestEffort(1);
  reanimator.RunBestEffort(2);

  // `Stop` waits for all the queued runs, including the best-effort ones, to
  // complete.
  reanimator.Stop();

  // Since the reanimator is started *before* queueing the actions, it's
  // possible that action 1 would finish before action 2 is queued.  Therefore,
  // we cannot assume a definite order here.
  EXPECT_THAT(processed, UnorderedElementsAre(1, 2));
}

// Checks that cancellation kills the right best-effort runs.
TEST_F(ReanimatorTest, CancelKillsBestEffortRuns) {
  std::vector<int> processed;
  absl::Notification proceed;
  ToyReanimator reanimator([&proceed, &processed](int const& key) {
    proceed.WaitForNotification();
    processed.push_back(key);
    return absl::OkStatus();
  });

  reanimator.Start();
  std::this_thread::sleep_for(100ms);

  // Queue two guaranteed runs that get stuck until `proceed` is notified.  At
  // some point, action 10 executes and action 5 is pending.
  reanimator.RunGuaranteed(10);
  reanimator.RunGuaranteed(5);

  // Queue best-effort runs.  They cannot proceed.
  reanimator.RunBestEffort(1);
  reanimator.RunBestEffort(2);
  reanimator.RunBestEffort(3);

  // Cancel the best-effort runs with a key strictly less than 3.
  reanimator.Cancel(/*before_key=*/3);

  // Unblock all the actions.
  proceed.Notify();

  reanimator.Stop();

  // The order is deterministic because all runs were queued and blocked at the
  // same time.
  EXPECT_THAT(processed, ElementsAre(10, 5, 3));
}

// Checks that a best-effort action that has a `RETURN_IF_STOPPED` observes the
// cancellation.
TEST_F(ReanimatorTest, CancelEverythingReturnIfStopped) {
  std::vector<int> processed;
  absl::Notification proceed1;
  absl::Notification proceed2;
  ToyReanimator reanimator([&proceed1, &proceed2, &processed](int const& key) {
    proceed1.WaitForNotification();
    RETURN_IF_STOPPED;
    proceed2.WaitForNotification();
    processed.push_back(key);
    return absl::OkStatus();
  });

  reanimator.Start();

  // Queue best-effort runs.  They cannot proceed.
  reanimator.RunBestEffort(1);
  reanimator.RunBestEffort(2);
  reanimator.RunBestEffort(3);

  // Wait for action 3 to actually run.
  std::this_thread::sleep_for(100ms);

  // This thread would cancel all the best-effort runs, but it cannot finish
  // while action 3 is running.
  std::thread cancellator([&reanimator]() {
    reanimator.Cancel(/*before_key=*/6);
  });

  // Wait for `Cancel` to get stuck.
  std::this_thread::sleep_for(100ms);

  // Unblock action 3.  It sees that it is stopped so it returns without ever
  // waiting for `proceed2`.
  proceed1.Notify();

  reanimator.Stop();

  // No action was run to completion.
  EXPECT_THAT(processed, IsEmpty());
}

TEST_F(ReanimatorTest, CancelEverythingSlow) {
  std::vector<int> processed;
  absl::Notification running;
  ToyReanimator reanimator([&processed, &running](int const& key) {
    running.Notify();
    std::this_thread::sleep_for(1000ms);
    processed.push_back(key);
    return absl::OkStatus();
  });

  // Queue best-effort runs.
  reanimator.RunBestEffort(1);
  reanimator.RunBestEffort(2);
  reanimator.RunBestEffort(3);

  // Start the reanimator and wait until action 3 is actually running.
  reanimator.Start();
  running.WaitForNotification();

  // Cancel all the best effort runs.  The action does not observe the
  // cancellation so `Cancel` only returns once it finishes.
  reanimator.Cancel(/*before_key=*/6);

  reanimator.Stop();

  // No action was run to completion.
  EXPECT_THAT(processed, ElementsAre(3));
}

TEST_F(ReanimatorTest, WaitWithProgressCallback) {
  absl::Mutex lock;
  std::vector<int> processed;
  ToyReanimator reanimator([&lock, &processed](int const& key) {
    absl::MutexLock l(&lock);
    processed.push_back(key);
    return absl::OkStatus();
  });
  reanimator.Start();

  auto const handle1 = reanimator.RunGuaranteed(1);
  auto const handle2 = reanimator.RunGuaranteed(2);

  std::vector<int> callback_keys;
  absl::Mutex callback_lock;
  EXPECT_OK(reanimator.Wait(
      handle2,
      [&callback_lock, &callback_keys](int const& key,
                                       absl::Status const& status) {
        absl::MutexLock l(&callback_lock);
        callback_keys.push_back(key);
      }));

  EXPECT_OK(reanimator.Wait(handle1));
  reanimator.Stop();

  // The progress callback was called at least for the run being waited for.
  absl::MutexLock l(&callback_lock);
  EXPECT_THAT(callback_keys, ElementsAre(2));
}

}  // namespace base
}  // namespace principia