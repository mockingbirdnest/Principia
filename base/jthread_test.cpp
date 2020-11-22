#include "base/jthread.hpp"

#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

TEST(JThreadTest, StopOnDestruction) {
  auto f = [](stop_token st, int value, bool* const observed_stop) {
    while (!st.stop_requested()) {
      std::cout << ++value << ' ' << std::flush;
      absl::SleepFor(absl::Milliseconds(10));
    }
    std::cout << std::endl;
    *observed_stop = true;
  };

  bool observed_stop = false;
  {
    jthread thread(f, 5, &observed_stop);
    absl::SleepFor(absl::Milliseconds(100));
  }
  EXPECT_TRUE(observed_stop);
}

TEST(JThreadTest, RequestStop) {
  bool observed_stop = false;

  jthread sleepy_worker(
      [](stop_token st, bool* const observed_stop) {
        for (int i = 0; i < 10; i++) {
          absl::SleepFor(absl::Milliseconds(10));
          if (st.stop_requested()) {
            std::cout << "Sleepy worker is requested to stop\n";
            *observed_stop = true;
            return;
          }
          std::cout << "Sleepy worker goes back to sleep\n";
        }
      },
      &observed_stop);

  absl::SleepFor(absl::Milliseconds(30));
  std::cout << "Requesting stop of sleepy worker\n";
  sleepy_worker.request_stop();
  sleepy_worker.join();
  EXPECT_TRUE(observed_stop);
}

TEST(JThreadTest, StopCallback) {
  // A worker thread.  It will wait until it is requested to stop.
  jthread worker([](stop_token st) {
    while (!st.stop_requested()) {
      std::cout << "Worker thread's id: " << std::this_thread::get_id() << '\n';
      absl::SleepFor(absl::Milliseconds(10));
    }
  });

  // Register a stop callback on the worker thread.
  bool callback_executed = false;
  stop_callback callback(worker.get_stop_token(), [&callback_executed] {
    std::cout << "Stop callback executed by thread: "
              << std::this_thread::get_id() << '\n';
    callback_executed = true;
  });

  // stop_callback objects can be destroyed prematurely to prevent execution.
  bool scoped_callback_executed = false;
  {
    stop_callback scoped_callback(
        worker.get_stop_token(), [&scoped_callback_executed] {
          // This will not be executed.
          std::cout << "Scoped stop callback executed by thread: "
                    << std::this_thread::get_id() << '\n';
          scoped_callback_executed = true;
        });
  }

  // Demonstrate which thread executes the stop_callback and when.  Define a
  // stopper function.
  int stop_count = 0;
  auto stopper_fn = [&stop_count, &worker] {
    if (worker.request_stop()) {
      std::cout << "Stop request executed by thread: "
                << std::this_thread::get_id() << '\n';
      ++stop_count;
    } else {
      std::cout << "Stop request not executed by thread: "
                << std::this_thread::get_id() << '\n';
    }
  };

  // Let multiple threads compete for stopping the worker thread.
  std::thread stopper1(stopper_fn);
  std::thread stopper2(stopper_fn);
  stopper1.join();
  stopper2.join();

  EXPECT_TRUE(callback_executed);
  EXPECT_FALSE(scoped_callback_executed);
  EXPECT_EQ(1, stop_count);

  // After a stop has already been requested, a new stop_callback executes
  // immediately.
  std::cout << "Main thread: " << std::this_thread::get_id() << '\n';

  bool callback_after_stop_executed = false;
  stop_callback callback_after_stop(
      worker.get_stop_token(), [&callback_after_stop_executed] {
        std::cout << "Stop callback executed by thread: "
                  << std::this_thread::get_id() << '\n';
        callback_after_stop_executed = true;
      });
  EXPECT_TRUE(callback_after_stop_executed);
}

TEST(JThreadTest, ThisJThread) {
  bool observed_stop = false;

  auto sleepy_worker = this_jthread::Make(
      [](stop_token /*unused*/, bool* const observed_stop) {
        for (int i = 0; i < 10; i++) {
          absl::SleepFor(absl::Milliseconds(10));
          if (this_jthread::get_stop_token().stop_requested()) {
            std::cout << "Sleepy worker is requested to stop\n";
            *observed_stop = true;
            return;
          }
          std::cout << "Sleepy worker goes back to sleep\n";
        }
      },
      &observed_stop);

  absl::SleepFor(absl::Milliseconds(30));
  std::cout << "Requesting stop of sleepy worker\n";
  sleepy_worker.request_stop();
  sleepy_worker.join();
  EXPECT_TRUE(observed_stop);
}



}  // namespace base
}  // namespace principia
