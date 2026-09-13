#include "base/reanimator.hpp"

#include <chrono>
#include <cstdint>
#include <string>
#include <thread>
#include <vector>

#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"
#include "absl/synchronization/notification.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace base {

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;
using namespace principia::base::_reanimator;
using namespace std::chrono_literals;

class ReanimatorTest : public ::testing::Test {
 protected:
  using ToyReanimator = Reanimator<int>;
};

TEST_F(ReanimatorTest, StartStop) {
  ToyReanimator reanimator([](int const) {
    return absl::OkStatus();
  });

  // May stop before starting.
  reanimator.Stop();

  // Idempotence.
  reanimator.Start();
  reanimator.Start();
  reanimator.Stop();
  reanimator.Stop();
}

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
  ToyReanimator reanimator([&processed](int const key) {
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

TEST_F(ReanimatorTest, CancelNoRun) {
  ToyReanimator reanimator([](int const) {
    return absl::OkStatus();
  });

  reanimator.Start();
  std::this_thread::sleep_for(100ms);
  reanimator.Cancel(/*before_key=*/1);
  reanimator.Stop();
}

// Checks that cancellation kills the right best-effort runs.
TEST_F(ReanimatorTest, CancelKillsBestEffortRuns) {
  std::vector<int> processed;
  absl::Notification proceed;
  ToyReanimator reanimator([&proceed, &processed](int const key) {
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
  ToyReanimator reanimator([&proceed1, &proceed2, &processed](int const key) {
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

  cancellator.join();
}

TEST_F(ReanimatorTest, CancelEverythingSlow) {
  std::vector<int> processed;
  absl::Notification running;
  ToyReanimator reanimator([&processed, &running](int const key) {
    running.Notify();
    std::this_thread::sleep_for(500ms);
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
  std::vector<int> processed;
  ToyReanimator reanimator([&processed](int const key) {
    processed.push_back(key);
    return absl::OkStatus();
  });

  // Queue three runs.
  auto const handle1 = reanimator.RunGuaranteed(1);
  reanimator.RunBestEffort(2);
  reanimator.RunGuaranteed(3);

  // Delay starting the reanimator until after the call to `Wait` to ensure
  // that the progress callback is called for all three actions.
  auto starter = std::thread([&reanimator]() {
    std::this_thread::sleep_for(100ms);
    reanimator.Start();
  });

  // Action 1 is the last executed, so the progress callback is called for all
  // three actions.
  std::vector<int> callback_keys;
  EXPECT_OK(reanimator.Wait(
      handle1, [&callback_keys](int const key, absl::Status const& status) {
        CHECK_OK(status);
        callback_keys.push_back(key);
      }));

  reanimator.Stop();

  EXPECT_THAT(callback_keys, ElementsAre(3, 2, 1));

  starter.join();
}

TEST_F(ReanimatorTest, Parameters) {
  std::vector<int> keys;
  std::vector<std::string> strings;
  std::vector<bool> bools;
  Reanimator<int, std::string, bool> reanimator(
      [&bools, &keys, &strings](
          int const key, std::string const& s, bool const b) {
        keys.push_back(key);
        strings.push_back(s);
        bools.push_back(b);
        return absl::OkStatus();
      });

  // Queue three runs.
  auto const handle1 = reanimator.RunGuaranteed(1, "1", true);
  reanimator.RunBestEffort(2, "2", false);
  auto const handle3 = reanimator.RunGuaranteed(3, "3", true);

  reanimator.Start();

  EXPECT_OK(reanimator.Wait(handle1));
  EXPECT_OK(reanimator.Wait(handle3));

  reanimator.Stop();

  EXPECT_THAT(keys, UnorderedElementsAre(3, 2, 1));
  EXPECT_THAT(strings, UnorderedElementsAre("3", "2", "1"));
  EXPECT_THAT(bools, UnorderedElementsAre(true, false, true));
}

}  // namespace base
}  // namespace principia
