export module principia.base.graveyard;

import std;

import principia.base.thread_pool;
import principia.base.graveyard;

namespace principia {
namespace base {
namespace _graveyard {
namespace internal {

using namespace principia::base::_thread_pool;

export {

class Graveyard {
 public:
  explicit Graveyard(std::int64_t number_of_threads);

  template<std::movable T>
  void Bury(T t);

 private:
  ThreadPool<void> gravedigger_;
};

}  // export
}  // namespace internal

export using internal::Graveyard;

}  // namespace _graveyard
}  // namespace base
}  // namespace principia

#include "base/graveyard_body.hpp"
