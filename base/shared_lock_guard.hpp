#pragma once

#include "base/macros.hpp"
#include "base/not_null.hpp"

#include <shared_mutex>

namespace principia {
namespace base {

// A helper class that the language designer didn't think useful of providing.
template<typename Mutex>
class shared_lock_guard final {
 public:
  explicit shared_lock_guard(Mutex& mutex);
  ~shared_lock_guard();

 private:
  not_null<Mutex*> const mutex_;
};

}  // namespace base
}  // namespace principia

#include "base/shared_lock_guard_body.hpp"
