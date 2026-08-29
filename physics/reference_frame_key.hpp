#pragma once

#include <cstddef>
#include <memory>
#include <utility>

#include "base/not_null.hpp"
#include "physics/ephemeris.hpp"
#include "physics/reference_frame.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _reference_frame_key {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_reference_frame;

// A copyable container for using `ReferenceFrame` and its subclasses as keys in
// hash maps.
template<typename InertialFrame, typename ThisFrame>
class ReferenceFrameKey {
 public:
  explicit ReferenceFrameKey(
      not_null<std::unique_ptr<ReferenceFrame<InertialFrame, ThisFrame> const>>
          frame);

  void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const;

  static ReferenceFrameKey ReadFromMessage(
      serialization::ReferenceFrame const& message,
      not_null<Ephemeris<InertialFrame> const*> ephemeris);

  template<typename H>
  friend H AbslHashValue(H h, ReferenceFrameKey const& key) {
    key.frame_->HashValue(absl::HashState::Create(&h));
    return std::move(h);
  }

  friend bool operator==(ReferenceFrameKey const& lhs,
                         ReferenceFrameKey const& rhs) {
    return *lhs.frame_ == *rhs.frame_;
  }

  // Support for heterogeneous lookups.
  struct absl_container_hash {
    using F = ReferenceFrame<InertialFrame, ThisFrame>;
    using is_transparent = void;
    std::size_t operator()(not_null<F const*> frame) const;
    std::size_t operator()(ReferenceFrameKey const& key) const;
  };

  struct absl_container_eq {
    using F = ReferenceFrame<InertialFrame, ThisFrame>;
    using is_transparent = void;
    bool operator()(ReferenceFrameKey const& lhs,
                    ReferenceFrameKey const& rhs) const;
    bool operator()(not_null<F const*> lhs, ReferenceFrameKey const& rhs) const;
    bool operator()(ReferenceFrameKey const& lhs, not_null<F const*> rhs) const;
  };

 private:
  not_null<std::shared_ptr<ReferenceFrame<InertialFrame, ThisFrame> const>>
      frame_;
};

}  // namespace internal

using internal::ReferenceFrameKey;

}  // namespace _reference_frame_key
}  // namespace physics
}  // namespace principia

#include "physics/reference_frame_key_body.hpp"
