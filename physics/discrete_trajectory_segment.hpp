#pragma once

#include <cstdint>

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory_iterator.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/discrete_trajectory_types.hpp"

namespace principia {
namespace physics {

class DiscreteTrajectoryIteratorTest;

namespace internal_discrete_trajectory_segment {

using geometry::Instant;
using physics::DegreesOfFreedom;

template<typename Frame>
class DiscreteTrajectorySegment {
 public:
  DiscreteTrajectorySegment() = default;
  virtual ~DiscreteTrajectorySegment() = default;

  // Moveable.
  DiscreteTrajectorySegment(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment& operator=(DiscreteTrajectorySegment&&) = default;
  DiscreteTrajectorySegment(const DiscreteTrajectorySegment&) = delete;
  DiscreteTrajectorySegment& operator=(const DiscreteTrajectorySegment&) =
      delete;

  DiscreteTrajectoryIterator<Frame> begin() const;
  DiscreteTrajectoryIterator<Frame> end() const;

  DiscreteTrajectoryIterator<Frame> rbegin() const;
  DiscreteTrajectoryIterator<Frame> rend() const;

  DiscreteTrajectoryIterator<Frame> find(Instant const& t) const;

  DiscreteTrajectoryIterator<Frame> lower_bound(Instant const& t) const;
  DiscreteTrajectoryIterator<Frame> upper_bound(Instant const& t) const;

  bool empty() const;
  virtual std::int64_t size() const;

 protected:
  // For mocking.
  using Timeline = internal_discrete_trajectory_types::Timeline<Frame>;

 private:
  void Append(Instant const& t,
              DegreesOfFreedom<Frame> const& degrees_of_freedom);

  void ForgetAfter(Instant const& t);
  void ForgetAfter(typename Timeline::const_iterator begin);

  void ForgetBefore(Instant const& t);
  void ForgetBefore(typename Timeline::const_iterator end);

  virtual typename Timeline::const_iterator timeline_begin() const;
  virtual typename Timeline::const_iterator timeline_end() const;

  DiscreteTrajectorySegmentIterator<Frame> self_;

  Timeline timeline_;
  absl::btree_set<Instant> dense_points_;

  template<typename F>
  friend class internal_discrete_trajectory_iterator::
      DiscreteTrajectoryIterator;
  friend class DiscreteTrajectoryIteratorTest;
};

}  // namespace internal_discrete_trajectory_segment

using internal_discrete_trajectory_segment::DiscreteTrajectorySegment;

}  // namespace physics
}  // namespace principia

#include "physics/discrete_trajectory_segment_body.hpp"
