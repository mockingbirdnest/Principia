#pragma once

#include "physics/reference_frame_key.hpp"

#include <memory>
#include <utility>

namespace principia {
namespace physics {
namespace _reference_frame_key {
namespace internal {

template<typename InertialFrame, typename ThisFrame>
ReferenceFrameKey<InertialFrame, ThisFrame>::ReferenceFrameKey(
    not_null<std::unique_ptr<ReferenceFrame<InertialFrame, ThisFrame>>> frame)
    : frame_(std::move(frame)) {}

template<typename H, typename InertialFrame, typename ThisFrame>
H AbslHashValue(H h, ReferenceFrameKey<InertialFrame, ThisFrame> const& key) {
  return H::combine(std::move(h), *key.frame_);
}

template<typename InertialFrame, typename ThisFrame>
bool operator==(ReferenceFrameKey<InertialFrame, ThisFrame> const& lhs,
                ReferenceFrameKey<InertialFrame, ThisFrame> const& rhs) {
  return *lhs.frame_ == *rhs.frame_;
}

}  // namespace internal
}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia
