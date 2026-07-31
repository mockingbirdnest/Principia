#pragma once

#include <memory>

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
  friend H AbslHashValue(H h, ReferenceFrameKey const& key);

 private:
  not_null<std::shared_ptr<ReferenceFrame<InertialFrame, ThisFrame>>> frame_;

  template<typename IF, typename TF>
  friend bool operator==(ReferenceFrameKey<IF, TF> const& lhs,
                         ReferenceFrameKey<IF, TF> const& rhs);
};

template<typename H, typename InertialFrame, typename ThisFrame>
H AbslHashValue(H h, ReferenceFrameKey<InertialFrame, ThisFrame> const& key);

template<typename InertialFrame, typename ThisFrame>
bool operator==(ReferenceFrameKey<InertialFrame, ThisFrame> const& lhs,
                ReferenceFrameKey<InertialFrame, ThisFrame> const& rhs);

}  // namespace internal

using internal::ReferenceFrameKey;

}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia

#include "physics/reference_frame_key_body.hpp"
