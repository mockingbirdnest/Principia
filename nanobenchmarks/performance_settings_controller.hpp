#pragma once

#include <memory>

#include "base/not_null.hpp"

namespace principia {
namespace nanobenchmarks {
namespace _performance_settings_controller {
namespace internal {

using namespace base::_not_null;

class PerformanceSettingsController {
 public:
  static not_null<std::unique_ptr<PerformanceSettingsController>> New();
  virtual ~PerformanceSettingsController() = default;
};

}  // namespace internal

using internal::PerformanceSettingsController;

}  // namespace _performance_settings_controller
}  // namespace nanobenchmarks
}  // namespace principia
