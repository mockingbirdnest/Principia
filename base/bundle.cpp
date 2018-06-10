
#include "base/bundle.hpp"

#include <functional>
#include <list>

#include "base/map_util.hpp"
#include "base/status.hpp"

#if !OS_MACOSX

namespace principia {
namespace base {

thread_local std::function<bool()> AbortRequested = [] { return false; };

Bundle::Bundle(int const workers)
    : max_workers_(workers),
      master_abort_(&AbortRequested) {
  workers_.reserve(max_workers_);
}

Status Bundle::Join() {
  {
    std::unique_lock<std::mutex> lock(lock_);
    CHECK(wait_on_empty_);
    wait_on_empty_.Flop();
  }
  tasks_not_empty_or_terminate_.notify_all();
  for (auto& thread : workers_) {
    thread.join();
  }
  std::shared_lock<std::shared_mutex> status_lock(status_lock_);
  return status_;
}

Status Bundle::JoinWithin(std::chrono::steady_clock::duration Δt) {
  return JoinBefore(std::chrono::steady_clock::now() + Δt);
}

Status Bundle::JoinBefore(std::chrono::steady_clock::time_point t) {
  {
    std::unique_lock<std::shared_mutex> status_lock(status_lock_);
    deadline_ = t;
  }
  return Join();
}

void Bundle::Add(Task task) {
  {
    std::unique_lock<std::mutex> lock(lock_);
    CHECK(wait_on_empty_);
    if (workers_.size() < max_workers_) {
      workers_.emplace_back(&Bundle::Toil, this);
    }
    tasks_.emplace(task);
  }
  tasks_not_empty_or_terminate_.notify_one();
}

void Bundle::Toil() {
  Task current_task;
  AbortRequested = std::bind(&Bundle::BundleShouldAbort, this);
  for (;;) {
    {
      std::unique_lock<std::mutex> lock(lock_);
      tasks_not_empty_or_terminate_.wait(
          lock,
          [this] {
            return !tasks_.empty() || !wait_on_empty_ || Aborting();
          });
      // The call to |BundleShouldAbort| checks for deadline expiry and master
      // abort.
      if (tasks_.empty() || BundleShouldAbort()) {
        return;
      }
      current_task = std::move(tasks_.front());
      tasks_.pop();
    }
    Status const status = current_task();
    if (!status.ok()) {
      Abort(status);
    }
  }
}

// On the |workers_|, |AbortRequested| is set to |BundleShouldAbort|.  That
// function checks whether the |Bundle| is already |Aborting()|, and
// additionally checks whether the |deadline_| of the |Bundle| has expired, if
// any, and whether |(*master_abort_)()|; if either of those holds, |Abort| is
// called with an appropriate |Status|.
// |*master_abort| is |AbortRequested| on the master thread, thus ensuring that
// a master abort trickles down.  In order for this to happen even if the
// |tasks_| do not call |AbortRequested|, |BundleShouldAbort| is also called
// before dequeuing a new task.
bool Bundle::BundleShouldAbort() {
  if (!Aborting()) {
    if ((*master_abort_)()) {
      Abort(Status(Error::ABORTED, "abort requested on bundle master"));
    } else if (DeadlineExceeded()) {
      Abort(Status(Error::DEADLINE_EXCEEDED, "bundle deadline exceeded"));
    }
  }
  return Aborting();
}

void Bundle::Abort(Status const status) {
  std::unique_lock<std::shared_mutex> status_lock(status_lock_);
  CHECK(!status.ok());
  if (!status_.ok()) {
    // Already aborting.
    return;
  }
  status_ = status;
  tasks_not_empty_or_terminate_.notify_all();
}

bool Bundle::Aborting() {
  std::shared_lock<std::shared_mutex> status_lock(status_lock_);
  return !status_.ok();
}

bool Bundle::DeadlineExceeded() {
  std::shared_lock<std::shared_mutex> status_lock(status_lock_);
  return deadline_ && std::chrono::steady_clock::now() > *deadline_;
}

}  // namespace base
}  // namespace principia


#endif
