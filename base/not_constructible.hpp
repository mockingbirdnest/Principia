#pragma once

namespace principia {
namespace base {

// Type traits should inherit from this, so that they may use inheritance
// without virtual destructors.
struct not_constructible {
  not_constructible() = delete;
  ~not_constructible() = delete;
};

}  // namespace base
}  // namespace principia
