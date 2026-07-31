#pragma once

#include "physics/reference_frame_key.hpp"

#include <cstddef>
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

template<typename InertialFrame, typename ThisFrame>
size_t
ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_hash::operator()(
    F const* const frame) const {
  return absl::HashOf(*frame);
}

template<typename InertialFrame, typename ThisFrame>
size_t
ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_hash::operator()(
    ReferenceFrameKey const& key) const {
  return absl::HashOf(*frame);
}

}  // namespace internal
}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia
