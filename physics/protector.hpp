
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

// The protector helps with preventing asynchronous changes to timelines, i.e.,
// classes that associate some data with distinct instants.  It makes it
// possible to protect time ranges of the form [t_min, +∞[.  Any client that
// would want to touch the time range ]-∞, t[ should do so through
// RunWhenUnprotected, and the change will be delayed until ]-∞, t[ becomes
// unprotected.  This class is thread-safe.
class Protector {
 public:
  // A callback that may be run immediately or in a delayed manner when the
  // state of the protector permits it.
  using Callback = std::function<void()>;

  // If the range ]-∞, t[ is unprotected, |callback| is run immediately and this
  // function returns true.  Otherwise |callback| is delayed and will be run as
  // soon as ]-∞, t[ becomes unprotected; returns false in this case.  The
  // callbacks are run without any lock held, in time order.
  bool RunWhenUnprotected(Instant const& t, Callback callback);

  // Protects and unprotects the time range [t_min, +∞[.
  void Protect(Instant const& t_min);
  void Unprotect(Instant const& t_min);

 private:
  absl::Mutex lock_;
  std::multimap<Instant, Callback> callbacks_ GUARDED_BY(lock_);
  std::multiset<Instant> protection_start_times_ GUARDED_BY(lock_);
};

}  // namespace internal_protector

using internal_protector::Protector;

}  // namespace physics
}  // namespace principia
