
#pragma once

namespace principia {
namespace base {

// Type traits should inherit from this, so that they may use inheritance
// without virtual destructors.
struct type_trait {
  ~type_trait() = delete;
};

}  // namespace base
}  // namespace principia
