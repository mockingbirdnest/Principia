#pragma once

#include <memory>
#include <utility>

#include "base/not_null.hpp"
#include "physics/reference_frame.hpp"

namespace principia {
namespace physics {
namespace _reference_frame_key {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::physics::_reference_frame;

// A copyable container for using `ReferenceFrame` and its subclasses as keys in
// hash maps.
template<typename InertialFrame, typename ThisFrame>
class ReferenceFrameKey {
 public:
  explicit ReferenceFrameKey(
      not_null<std::unique_ptr<ReferenceFrame<InertialFrame, ThisFrame>>>
          frame);

  template<typename H>
  friend H AbslHashValue(H h, ReferenceFrameKey const& key) {
    return H::combine(std::move(h), *key.frame_);
  }

  friend bool operator==(ReferenceFrameKey const& lhs,
                         ReferenceFrameKey const& rhs) {
    return *lhs.frame_ == *rhs.frame_;
  }

 private:
  not_null<std::shared_ptr<ReferenceFrame<InertialFrame, ThisFrame>>> frame_;
};

}  // namespace internal

using internal::ReferenceFrameKey;

}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia

#include "physics/reference_frame_key_body.hpp"
