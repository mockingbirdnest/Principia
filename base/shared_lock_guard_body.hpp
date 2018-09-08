#pragma once

#include "base/shared_lock_guard.hpp"

namespace principia {
namespace base {

template<typename Mutex>
shared_lock_guard<Mutex>::shared_lock_guard(Mutex& mutex) : mutex_(&mutex) {
  mutex_->lock_shared();
}

template<typename Mutex>
shared_lock_guard<Mutex>::~shared_lock_guard() {
  mutex_->unlock_shared();
}

}  // namespace base
}  // namespace principia
