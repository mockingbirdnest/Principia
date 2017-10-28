#pragma once

#include "base/shared_lock_guard.hpp"

namespace principia {
namespace base {

#if HAS_SHARED_MUTEX

template<typename Mutex>
shared_lock_guard<Mutex>::shared_lock_guard(Mutex& mutex) : mutex_(&mutex) {
  mutex_->lock_shared();
}

template<typename Mutex>
shared_lock_guard<Mutex>::~shared_lock_guard() {
  mutex_->unlock_shared();
}

#endif

}  // namespace base
}  // namespace principia
