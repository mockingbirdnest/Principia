
#pragma once

#include <cstdint>
#include <memory>

#include "base/thread_pool.hpp"

namespace principia {
namespace base {

template<typename T>
class Graveyard {
 public:
  explicit Graveyard(std::int64_t number_of_threads);

  void Bury(std::unique_ptr<T> t);

 private:
  ThreadPool<void> gravedigger_;
};

}  // namespace base
}  // namespace principia

#include "base/graveyard_body.hpp"
