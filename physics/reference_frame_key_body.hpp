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
    not_null<std::unique_ptr<ReferenceFrame<InertialFrame, ThisFrame> const>>
        frame)
    : frame_(std::move(frame)) {}

template<typename InertialFrame, typename ThisFrame>
std::size_t
ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_hash::operator()(
    not_null<F const*> const frame) const {
  return absl::HashOf(*frame);
}

template<typename InertialFrame, typename ThisFrame>
std::size_t
ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_hash::operator()(
    ReferenceFrameKey const& key) const {
  return absl::HashOf(key);
}

template<typename InertialFrame, typename ThisFrame>
bool ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_eq::
operator()(ReferenceFrameKey const& lhs,
           ReferenceFrameKey const& rhs) const {
  return *lhs.frame_ == *rhs.frame_;
}

template<typename InertialFrame, typename ThisFrame>
bool ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_eq::
operator()(not_null<F const*> const lhs, ReferenceFrameKey const& rhs) const {
  return *lhs == *rhs.frame_;
}

template<typename InertialFrame, typename ThisFrame>
bool ReferenceFrameKey<InertialFrame, ThisFrame>::absl_container_eq::
operator()(ReferenceFrameKey const& lhs, not_null<F const*> const rhs) const {
  return *lhs.frame_ == *rhs;
}

}  // namespace internal
}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia
