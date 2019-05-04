
#pragma once

#include <functional>
#include <set>
#include <map>

#include "absl/base/thread_annotations.h"
#include "absl/synchronization/mutex.h"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_protector {

using geometry::Instant;

//TODO(phl):comments
// Holds the state of all the guards and the callbacks whose execution has been
// delayed due to protection by one or several guards.  This class is
// thread-safe.
class Protector {
 public:
  // A callback that may be run immediately or in a delayed manner when the
  // state of the guards permits it.
  using Callback = std::function<void()>;

  // If the range ]-∞, t[ is unprotected, |callback| is run immediately and this
  // function returns true.  Otherwise |callback| is delayed and will be run as
  // soon as ]-∞, t[ becomes unprotected; returns false in this case.  The
  // callback are run without any lock held, in time order.
  bool RunWhenUnprotected(Instant const& t, Callback callback);

  // Protects and unprotects the time range [t_min, +∞[.
  void Protect(Instant const& t_min);
  void Unprotect(Instant const& t_min);

 private:
  absl::Mutex lock_;
  std::multimap<Instant, Callback> callbacks_ GUARDED_BY(lock_);
  std::multiset<Instant> guard_start_times_ GUARDED_BY(lock_);
};

}  // namespace internal_protector

using internal_protector::Protector;

}  // namespace physics
}  // namespace principia
