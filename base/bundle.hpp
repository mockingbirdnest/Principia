#pragma once

#include <atomic>
#include <experimental/optional>
#include <functional>
#include <list>
#include <memory>
#include <mutex>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

thread_local std::atomic<bool> AbortRequested;

class Bundle {
  class Task;

 public:
  Bundle(int threads);

  // Returns the first non-OK status encountered.  All threads are joined;
  // no calls to |Add| or |Join| may follow this call.
  Status Join();

  // The |Task| returned by |Add| may be used in calls to |Abort| or |JoinTask|;
  // |this| shares ownership until |task| has returned.
  std::shared_ptr<Task> Add(std::function<Status()> task);

  // Removes the given |task| from the queue if it has not started, or instructs
  // it to cooperatively abort if it is running.
  void Abort(not_null<Task*> task);

  // Waits until the given |task| is scheduled and completes, and returns its
  // status.
  Status JoinTask(not_null<Task*> task);

 private:
  class Task {
   public:
    Task(std::function<Status()> body);

    std::function<Status()> body_;

    // |status_|, |in_queue_|, and |abort_requested_| are guarded by the |lock_|
    // of the |Bundle|.

    // Has a value if and only if the task has completed.
    std::experimental::optional<Status> status_;
    // Has a value if and only if the task has yet to start.
    std::experimental::optional<std::list<std::shared_ptr<Task>>::iterator>
        in_queue_;
    // Is non-null if and only if the task is running.
    std::atomic<bool>* abort_requested_ = nullptr;

    std::condition_variable done_;
  };

  std::mutex lock_;
  std::condition_variable tasks_not_empty_or_joining_;

  bool joining_ GUARDED_BY(lock_);
  std::list<std::shared_ptr<Task>> tasks_ GUARDED_BY(lock_);
  Status status_ GUARDED_BY(lock_);

  std::vector<std::thread> threads_;
};

}  // namespace base
}  // namespace principia
