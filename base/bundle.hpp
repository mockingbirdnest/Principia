#pragma once

#include <experimental/optional>
#include <functional>
#include <list>
#include <memory>
#include <mutex>
#include <vector>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"

namespace principia {
namespace base {

// Functions that can be cooperatively aborted should check |AbortRequested()|
// when appropriate.
extern thread_local std::function<bool()> AbortRequested;

class Bundle {
  class Task;
 public:
  class TaskHandle {
   public:
    TaskHandle() = default;
    TaskHandle(TaskHandle const&) = delete;
    TaskHandle(TaskHandle&&) = default;
    TaskHandle& operator=(TaskHandle const&) = delete;
    TaskHandle& operator=(TaskHandle&&) = default;

   private:
    TaskHandle(std::shared_ptr<Task> task);
    std::shared_ptr<Task> task_;
    friend class Bundle;
  };

  Bundle(int workers);

  // Returns the first non-OK status encountered.  All worker threads are
  // joined; no calls to |Add| or |Join| may follow this call.
  Status Join();

  // The |TaskHandle| returned by |Add| may be used in a call to |Abort| or
  // |JoinTask|.  |this| shares ownership of the |Task| with the handle until
  // the |Task| is |Done|.
  TaskHandle Add(std::function<Status()> task_body);

  // Removes the given |task| from the queue if it has not started, or instructs
  // it to cooperatively abort if it is running.
  void Abort(TaskHandle const& task);

  // Waits until the given |task| completes, and returns its status.  Takes
  // ownership of |task|, so that |*task->task_| is destroyed when this function
  // returns.
  Status JoinTask(TaskHandle const task);

 private:
  class Task {
   public:
    // All state transitions set additional bits; No state transition clears a
    // bit.
    enum ExecutionState : std::uint8_t {
      // |this| is in the queue waiting for a thread to execute it.
      Inactive = 0b000,
      // |body_()| is being executed.
      Active = 0b001,
      // |body_()| has returned or |this| has been dequeued from the |Inactive|
      // state by an |Abort|.
      Terminated = 0b010,
      // |Abort| has been called.
      Aborted = 0b100,
    };

    Task(not_null<Bundle*> bundle, std::function<Status()> body);

    std::uint8_t execution_state() const;

    // |bundle_->lock_| must be held by the caller.  Requires
    // |execution_state() == Inactive|.
    void set_queue_it(std::list<std::shared_ptr<Task>>::iterator it);

    // No prerequisites.  If |Inactive|, |body_()| never runs.  If |Active|,
    // a cooperative abort is requested.  Sets the |Aborted| bit.
    void Abort();

    // Requires |execution_state() == Inactive|.  Sets the |Active| bit, and
    // calls |body_()|.  When |body_()| returns, sets the |Terminated| bit, and
    // notifies |terminated_|.
    void Execute();

    // Requires |execution_state() & Terminated|.
    // Returns the status returned by |body()|, or a status with
    // |Error::ABORTED| if |Abort| was called while |Inactive|.
    Status status() const;

    // Returns when the |Terminated| bit is set.
    void wait_for_termination();

   private:
    not_null<Bundle*> const bundle_;
    std::function<Status()> body_;

    std::uint8_t execution_state_ GUARDED_BY(bundle_->lock_) = Inactive;
    Status status_ GUARDED_BY(bundle_->lock_);
    // If |execution_state_ == Inactive|, the position of this task in the
    // queue, used to remove it from the queue.
    std::list<std::shared_ptr<Task>>::iterator queue_it
        GUARDED_BY(bundle_->lock_);

    std::condition_variable terminated_;
  };

  // The |workers_| |Toil|.
  void Toil();

  std::mutex lock_;
  // Notified once when a task is available for execution.  Notified to all
  // waiting workers when |terminate_on_empty_| is set.  Workers wait on this
  // when no |tasks_| are available.
  std::condition_variable tasks_not_empty_or_terminate_;

  // Whether the workers should terminate when no tasks are available.  Set by
  // |Join|.
  bool terminate_on_empty_ GUARDED_BY(lock_);
  // The queue of |Tasks| (it is a |list| so that aborted tasks may be removed
  // from the queue).  Insertion at the |back|, extraction at the |front|.
  std::list<std::shared_ptr<Task>> tasks_ GUARDED_BY(lock_);
  // The first erroneous status of a |Terminated| task, or |OK| if everything is
  // fine.
  Status status_ GUARDED_BY(lock_);

  // These variables are only accessed by the thread that constructs the
  // |Bundle| and calls |Add|.
  std::vector<std::thread> workers_;
  int const max_workers_;
};

}  // namespace base
}  // namespace principia
