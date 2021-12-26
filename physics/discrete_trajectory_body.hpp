#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_cat.h"
#include "base/status_utilities.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory {

using base::make_not_null_unique;
using base::uninitialized;
using geometry::InfiniteFuture;
using geometry::InfinitePast;
using quantities::Length;

template<typename Frame>
DiscreteTrajectory<Frame>::DiscreteTrajectory()
    : segments_(make_not_null_unique<Segments>(1)) {
  auto const sit = segments_->begin();
  auto const self = SegmentIterator(segments_.get(), sit);
  *sit = DiscreteTrajectorySegment<Frame>(self);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::reference
DiscreteTrajectory<Frame>::front() const {
  return *begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::reference
DiscreteTrajectory<Frame>::back() const {
  return *rbegin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::iterator
DiscreteTrajectory<Frame>::begin() const {
  return segments_->front().begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::iterator
DiscreteTrajectory<Frame>::end() const {
  return segments_->back().end();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::reverse_iterator
DiscreteTrajectory<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::reverse_iterator
DiscreteTrajectory<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::empty() const {
  return segment_by_left_endpoint_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::size() const {
  if (empty()) {
    return 0;
  }
  std::int64_t size = 1;
  std::int64_t nonempty_segments = 0;
  for (auto const& segment : *segments_) {
    if (!segment.empty()) {
      ++nonempty_segments;
      size += segment.size();
    }
  }
  size -= nonempty_segments;  // The junction points.
  return size;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::clear() {
  segments_->erase(std::next(segments_->begin()), segments_->end());
  segments_->front().clear();
  segment_by_left_endpoint_.clear();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::iterator
DiscreteTrajectory<Frame>::find(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    return end();
  }
  auto const sit = leit->second;
  auto const it = sit->find(t);
  if (it == sit->end()) {
    return end();
  }
  return it;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::iterator
DiscreteTrajectory<Frame>::lower_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    // This includes an empty trajectory.
    return begin();
  }
  auto const sit = leit->second;
  auto const it = sit->lower_bound(t);
  if (it == sit->end()) {
    return end();
  }
  return it;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::iterator
DiscreteTrajectory<Frame>::upper_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    // This includes an empty trajectory.
    return begin();
  }
  auto const sit = leit->second;
  auto const it = sit->upper_bound(t);
  if (it == sit->end()) {
    return end();
  }
  return it;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::SegmentRange
DiscreteTrajectory<Frame>::segments() const {
  return SegmentRange(SegmentIterator(
                          segments_.get(), segments_->begin()),
                      SegmentIterator(
                          segments_.get(), segments_->end()));
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::ReverseSegmentRange
DiscreteTrajectory<Frame>::rsegments() const {
  return ReverseSegmentRange(std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->end())),
                             std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->begin())));
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::SegmentIterator
DiscreteTrajectory<Frame>::NewSegment() {
  auto& last_segment = segments_->back();

  segments_->emplace_back();
  auto const new_segment_sit = --segments_->end();
  auto const new_self = SegmentIterator(segments_.get(), new_segment_sit);
  *new_segment_sit = DiscreteTrajectorySegment<Frame>(new_self);

  // It is only possible to insert a segment after an empty segment if the
  // entire trajectory is empty.
  if (last_segment.empty()) {
    CHECK(segment_by_left_endpoint_.empty())
        << "Inserting after an empty segment but the trajectory is not empty, "
        << "its time-to-segment map has " << segment_by_left_endpoint_.size()
        << " entries";
  } else {
    // Duplicate the last point of the previous segment.
    auto const& [last_time, last_degrees_of_freedom] = *last_segment.rbegin();
    new_segment_sit->Append(last_time, last_degrees_of_freedom).IgnoreError();
    // The use of |insert_or_assign| ensure that we override any entry with the
    // same left endpoint.
    segment_by_left_endpoint_.insert_or_assign(
        segment_by_left_endpoint_.end(), last_time, new_segment_sit);
  }

  CHECK_OK(ConsistencyStatus());
  return new_self;
}

template<typename Frame>
DiscreteTrajectory<Frame>
DiscreteTrajectory<Frame>::DetachSegments(SegmentIterator const begin) {
  DiscreteTrajectory detached(uninitialized);

  // Move the detached segments to the new trajectory.
  detached.segments_->splice(detached.segments_->end(),
                             *segments_,
                             begin.iterator(), segments_->end());

  AdjustAfterSplicing(/*from=*/*this,
                      /*to=*/detached,
                      /*to_segments_begin=*/detached.segments_->begin());

  CHECK_OK(ConsistencyStatus());
  return detached;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::SegmentIterator
DiscreteTrajectory<Frame>::AttachSegments(DiscreteTrajectory trajectory) {
  CHECK(!trajectory.empty());

  if (empty()) {
    *this = DiscreteTrajectory(uninitialized);
  } else if (back().time == trajectory.front().time) {
    CHECK_EQ(back().degrees_of_freedom, trajectory.front().degrees_of_freedom)
        << "Mismatching degrees of freedom when attaching segments";
  } else {
    // If the points are not matching, prepend a matching point to |trajectory|
    // and update the time-to-segment map.
    CHECK_LT(back().time, trajectory.front().time)
        << "Mismatching times when attaching segments";
    trajectory.segments_->begin()->Prepend(back().time,
                                           back().degrees_of_freedom);
    auto const leit = trajectory.segment_by_left_endpoint_.cbegin();
    auto const sit = leit->second;
    trajectory.segment_by_left_endpoint_.erase(leit);
    trajectory.segment_by_left_endpoint_.emplace(back().time, sit);
  }

  // The |end| iterator keeps pointing at the end after the splice.  Instead,
  // we track the iterator to the last segment.
  auto const last_before_splice = --segments_->end();

  // Move the attached segments to the end of this trajectory.
  segments_->splice(segments_->end(),
                    *trajectory.segments_);

  auto const end_before_splice = std::next(last_before_splice);
  AdjustAfterSplicing(/*from=*/trajectory,
                      /*to=*/*this,
                      /*to_segments_begin=*/end_before_splice);

  CHECK_OK(ConsistencyStatus());
  return SegmentIterator(segments_.get(), end_before_splice);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::DeleteSegments(SegmentIterator& begin) {
  segments_->erase(begin.iterator(), segments_->end());
  if (segments_->empty()) {
    segment_by_left_endpoint_.clear();
  } else {
    auto const last_segment = --segments_->end();
    if (last_segment->empty()) {
      segment_by_left_endpoint_.clear();
    } else {
      // If there remains a non-empty segment, update the time-to-segment map
      // with for its time, and delete the rest.
      auto const& [leit, _] = segment_by_left_endpoint_.insert_or_assign(
          last_segment->front().time, last_segment);
      segment_by_left_endpoint_.erase(std::next(leit),
                                      segment_by_left_endpoint_.end());
    }
  }
  // Make sure that the client doesn't try to use the invalid iterator.
  begin = segments().end();

  CHECK_OK(ConsistencyStatus());
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& t) {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.end()) {
    clear();
    return;
  }
  auto const sit = leit->second;
  sit->ForgetAfter(t);
  // Here |sit| designates a segment starting at or after |t|.  If |t| is
  // exactly at the beginning of the segment,
  // |DiscreteTrajectorySegment::ForgetAfter| will leave it empty.  In that
  // case we drop the segment entirely, unless it is the only one in the
  // trajectory.
  if (sit->empty()) {
    if (sit == segments_->begin()) {
      segments_->erase(std::next(sit), segments_->end());
    } else {
      segments_->erase(sit, segments_->end());
    }
    segment_by_left_endpoint_.erase(leit, segment_by_left_endpoint_.end());
  } else {
    segments_->erase(std::next(sit), segments_->end());
    segment_by_left_endpoint_.erase(std::next(leit),
                                    segment_by_left_endpoint_.end());
  }

  CHECK_OK(ConsistencyStatus());
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(iterator const it) {
  if (it != end()) {
    ForgetAfter(it->time);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& t) {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.end()) {
    return;
  }
  auto const sit = leit->second;
  // This call may make the segment |*sit| empty if |t| is after the end of
  // |*sit|.
  // NOTE(phl): This declaration is necessary because MSVC corrupts |t| during
  // the call below.
  Instant const t_saved = t;
  sit->ForgetBefore(t);
  for (auto s = segments_->begin(); s != sit; ++s) {
    // This call may either make the segment |*s| empty or leave it with a
    // single point matching |sit->front()|.
    s->ForgetBefore(t_saved);
  }

  // Erase all the entries before and including |leit|.  These are entries for
  // now-empty segments or for 1-point segments, and the segment which we may
  // just have truncated.
  segment_by_left_endpoint_.erase(segment_by_left_endpoint_.begin(),
                                  std::next(leit));
  // It |*sit| is not empty, recreate an entry with its new left endpoint.
  if (!sit->empty()) {
    segment_by_left_endpoint_.insert_or_assign(
        segment_by_left_endpoint_.begin(),
        sit->front().time,
        sit);
  }

  CHECK_OK(ConsistencyStatus());
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(iterator const it) {
  if (it == end()) {
    clear();
  } else {
    ForgetBefore(it->time);
  }
}

template<typename Frame>
absl::Status DiscreteTrajectory<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  typename Segments::iterator sit;
  if (segment_by_left_endpoint_.empty()) {
    // If this is the first point appended to this trajectory, insert it in the
    // time-to-segment map.
    sit = --segments_->end();
    segment_by_left_endpoint_.insert_or_assign(
        segment_by_left_endpoint_.end(), t, sit);
  } else {
    auto const leit = FindSegment(t);
    CHECK(leit != segment_by_left_endpoint_.end())
        << "Append at " << t << " before the beginning of the trajectory at "
        << front().time;
    // The segment is expected to always have a point copied from its
    // predecessor.
    sit = leit->second;
    CHECK(!sit->empty()) << "Empty segment at " << t;
  }
  RETURN_IF_ERROR(sit->Append(t, degrees_of_freedom));

  DCHECK_OK(ConsistencyStatus());
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Merge(DiscreteTrajectory<Frame> trajectory) {
  auto sit_s = trajectory.segments_->begin();  // Source iterator.
  auto sit_t = segments_->begin();  // Target iterator.
  for (;;) {
    if (sit_s != trajectory.segments_->end() && sit_t != segments_->end()) {
      // Record the existing left endpoint to update the time-to-segment map as
      // needed.
      const std::optional<Instant> left_endpoint =
          sit_t->empty() ? std::nullopt
                         : std::make_optional(sit_t->front().time);

      // Merge corresponding segments.
      sit_t->Merge(std::move(*sit_s));

      // If the left endpoint of |sit_t| has changed, remove its entry from the
      // time-to-segment map, if any.
      if (left_endpoint.has_value() &&
          sit_t->front().time < left_endpoint.value()) {
        auto const it = segment_by_left_endpoint_.find(left_endpoint.value());
        if (it != segment_by_left_endpoint_.end() && it->second == sit_t) {
          segment_by_left_endpoint_.erase(left_endpoint.value());
        }
      }
      // Insert a new entry in the time-to-segment map if the segment is not
      // empty.  This entry will be overwritten by any future entry at the same
      // time, thereby enforcing the invariants of the time-to-segment map.
      if (!sit_t->empty()) {
        segment_by_left_endpoint_.insert_or_assign(sit_t->front().time, sit_t);
      }

      ++sit_s;
      ++sit_t;
    } else if (sit_s != trajectory.segments_->end()) {
      // No more segments in the target.  We splice the segments of the source.

      // The |end| iterator keeps pointing at the end after the splice.
      // Instead, we track the iterator to the last segment.
      auto const last_before_splice = --segments_->end();

      segments_->splice(segments_->end(),
                        *trajectory.segments_,
                        sit_s,
                        trajectory.segments_->end());

      auto const end_before_splice = std::next(last_before_splice);
      AdjustAfterSplicing(/*from=*/trajectory,
                          /*to=*/*this,
                          /*to_segments_begin=*/end_before_splice);
      break;
    } else {
      // No more segments in the source, or both lists done.
      break;
    }
  }

  CHECK_OK(ConsistencyStatus());
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_min() const {
  if (empty()) {
    return InfiniteFuture;
  }
  for (auto sit = segments_->begin();; ++sit) {
    if (!sit->empty()) {
      return sit->t_min();
    }
  }
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_max() const {
  if (empty()) {
    return InfinitePast;
  }
  return segments_->back().t_max();
}

template<typename Frame>
Position<Frame> DiscreteTrajectory<Frame>::EvaluatePosition(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {
  WriteToMessage(message, begin(), end(), tracked, exact);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    iterator const begin,
    iterator const end,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {
  // Construct a map to efficiently find if a segment must be tracked.  The
  // keys are pointers to segments in |tracked|, the values are the
  // corresponding indices.
  absl::flat_hash_map<DiscreteTrajectorySegment<Frame> const*, int>
      segment_to_position;
  for (int i = 0; i < tracked.size(); ++i) {
    if (tracked[i] != segments().end()) {
      segment_to_position.emplace(&*tracked[i], i);
    }
  }

  // Initialize the tracked positions to be able to recognize if some are
  // missing.
  message->mutable_tracked_position()->Resize(
      tracked.size(),
      serialization::DiscreteTrajectory::MISSING_TRACKED_POSITION);

  // Convert to instants the iterators that define the range to write.  This is
  // necessary to do lookups in each segment to obtain segment-specific
  // iterators.
  Instant const begin_time = begin == this->end() ? InfiniteFuture
                                                  : begin->time;
  Instant const end_time = end == this->end() ? InfiniteFuture
                                              : end->time;

  // The set of segments that intersect the range to write.
  absl::flat_hash_set<
      DiscreteTrajectorySegment<Frame> const*> intersecting_segments;
  bool intersect_range = false;

  // The position of a segment in the repeated field |segment|.
  int segment_position = 0;
  for (auto sit = segments_->begin();
       sit != segments_->end();
       ++sit, ++segment_position) {
    // Look up in |*sit| the instants that define the range to write.
    // |lower_bound| and |upper_bound| return the past-the-end-of-segment
    // iterator if no point exists after the given time, i.e., for the segments
    // that precede the intersection.
    auto const begin_time_it = sit->lower_bound(begin_time);
    auto const end_time_it = sit->lower_bound(end_time);

    // If the above iterators determine an empty range, and we have already seen
    // a segment that intersects the range to write, we are past the
    // intersection.  Skip all the remaining segments.
    if (intersect_range && begin_time_it == end_time_it) {
      break;
    }

    // If |*sit| contains a point at or after |begin_time|, it intersects the
    // range to write.
    intersect_range = begin_time_it != sit->end();
    if (intersect_range) {
      intersecting_segments.insert(&*sit);
    }

    // Note that we execute this call for the segments that precede the
    // intersection in order to write the correct structure of (empty) segments.
    sit->WriteToMessage(
        message->add_segment(), begin_time_it, end_time_it, exact);

    if (auto const position_it = segment_to_position.find(&*sit);
        position_it != segment_to_position.end()) {
      // The field |tracked_position| is indexed by the indices in |tracked|.
      // Its value is the position of a tracked segment in the field |segment|.
      message->set_tracked_position(position_it->second, segment_position);
    }
  }

  // Write the left endpoints by scanning them in parallel with the segments.
  int i = 0;
  auto sit1 = segments_->begin();
  for (auto const& [t, sit2] : segment_by_left_endpoint_) {
    while (sit1 != sit2) {
      ++sit1;
      ++i;
    }
    // Skip the segments that don't intersect the range to write, as we pretend
    // that they are empty.  Adjust the left endpoint to account for the segment
    // that may have been truncated on the left.
    if (intersecting_segments.contains(&*sit2)) {
      auto* const segment_by_left_endpoint =
          message->add_segment_by_left_endpoint();
      Instant const left_endpoint = std::max(t, begin_time);
      left_endpoint.WriteToMessage(
          segment_by_left_endpoint->mutable_left_endpoint());
      segment_by_left_endpoint->set_segment(i);
    }
  }
}

template<typename Frame>
template<typename F, typename>
DiscreteTrajectory<Frame>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked) {
  DiscreteTrajectory trajectory(uninitialized);

  bool const is_pre_hamilton = message.segment_size() == 0;
  if (is_pre_hamilton) {
    LOG_IF(WARNING, is_pre_hamilton)
        << "Reading pre-Hamilton DiscreteTrajectory";
    ReadFromPreHamiltonMessage(
        message, tracked, /*fork_point=*/std::nullopt, trajectory);
    CHECK_OK(trajectory.ConsistencyStatus());
    return trajectory;
  }

  // First restore the segments themselves.  |segment_iterators| will be used to
  // restore the tracked segments.
  std::vector<SegmentIterator> segment_iterators;
  segment_iterators.reserve(message.segment_size());
  for (auto const& serialized_segment : message.segment()) {
    trajectory.segments_->emplace_back();
    auto const sit = --trajectory.segments_->end();
    auto const self = SegmentIterator(trajectory.segments_.get(), sit);
    *sit = DiscreteTrajectorySegment<Frame>::ReadFromMessage(serialized_segment,
                                                             self);
    segment_iterators.push_back(self);
  }

  // Restore the tracked segments.
  CHECK_EQ(tracked.size(), message.tracked_position_size());
  for (int i = 0; i < message.tracked_position_size(); ++i) {
    int const tracked_position = message.tracked_position(i);
    if (tracked_position ==
        serialization::DiscreteTrajectory::MISSING_TRACKED_POSITION) {
      *tracked[i] = trajectory.segments().end();
    } else {
      *tracked[i] = segment_iterators[tracked_position];
    }
  }

  // Finally restore the left endpoints.
  int i = 0;
  auto sit = trajectory.segments_->begin();
  for (auto const& segment_by_left_endpoint :
       message.segment_by_left_endpoint()) {
    auto const t =
        Instant::ReadFromMessage(segment_by_left_endpoint.left_endpoint());
    while (segment_by_left_endpoint.segment() != i) {
      ++sit;
      ++i;
    }
    trajectory.segment_by_left_endpoint_.insert_or_assign(
        trajectory.segment_by_left_endpoint_.end(), t, sit);
  }

  CHECK_OK(trajectory.ConsistencyStatus());
  return trajectory;
}

template<typename Frame>
DiscreteTrajectory<Frame>::DiscreteTrajectory(uninitialized_t)
    : segments_(make_not_null_unique<Segments>()) {}

template<typename Frame>
typename DiscreteTrajectory<Frame>::SegmentByLeftEndpoint::iterator
DiscreteTrajectory<Frame>::FindSegment(
    Instant const& t) {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  if (it == segment_by_left_endpoint_.begin()) {
    // This includes an empty trajectory.
    return segment_by_left_endpoint_.end();
  } else {
    return --it;
  }
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::SegmentByLeftEndpoint::const_iterator
DiscreteTrajectory<Frame>::FindSegment(
    Instant const& t) const {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  if (it == segment_by_left_endpoint_.begin()) {
    // This includes an empty trajectory.
    return segment_by_left_endpoint_.cend();
  } else {
    return --it;
  }
}

template<typename Frame>
absl::Status DiscreteTrajectory<Frame>::ConsistencyStatus() const {
  if (segments_->size() < segment_by_left_endpoint_.size()) {
    return absl::InternalError(absl::StrCat("Size mismatch ",
                                            segments_->size(),
                                            " and ",
                                            segment_by_left_endpoint_.size()));
  }
  if (segment_by_left_endpoint_.empty()) {
    int i = 0;
    for (auto const& segment : *segments_) {
      if (!segment.empty()) {
        return absl::InternalError(absl::StrCat(
            "Segment #", i,
            " is not empty but is not in the time-to-segment map"));
      }
      ++i;
    }
  }
  {
    int i = 0;
    for (auto const& [left_endpoint, sit] : segment_by_left_endpoint_) {
      if (sit->empty()) {
      } else if (left_endpoint != sit->front().time) {
        return absl::InternalError(
            absl::StrCat("Times mismatch ",
                          DebugString(left_endpoint),
                          " and ",
                          DebugString(sit->front().time),
                          " between segment #", i,
                          " and the time-to-segment map"));
      }
    }
  }
  {
    int i = 0;
    auto sit1 = segments_->cbegin();
    for (auto leit = segment_by_left_endpoint_.cbegin();
         leit != segment_by_left_endpoint_.cend();) {
      auto const sit2 = leit->second;
      if (sit2->empty()) {
        return absl::InternalError(
            absl::StrCat("Empty segment for time-to-segment entry #", i));
      } else if (sit1 == sit2) {
        ++sit1;
        ++leit;
        ++i;
      } else if (sit1->empty() || sit1->size() == 1) {
        ++sit1;
      } else {
        return absl::InternalError(
            absl::StrCat("Mismatch for time-to-segment entry #", i,
                         " at times ", DebugString(sit1->front().time),
                         " and ", DebugString(sit2->front().time)));
      }
    }
    if (sit1 != segments_->cend() && !sit1->empty()) {
        return absl::InternalError(
            absl::StrCat("Segment at time ", DebugString(sit1->front().time),
                         " missing in the time-to-segment map of size ",
                         segment_by_left_endpoint_.size()));
    }
  }
  if (!segments_->empty()) {
    int i = 0;
    for (auto sit = segments_->cbegin();
         sit != std::prev(segments_->cend());
         ++sit, ++i) {
      // Great care is required here because the DiscreteTrajectoryIterator will
      // "helpfully" paper over differences in degrees of freedom as long as the
      // times match.  We must look at the endpoints of the timeline explicitly.
      if (!sit->timeline_empty()) {
        auto const timeline_rbegin = --sit->timeline_end();
        auto const timeline_begin = std::next(sit)->timeline_begin();
        if (timeline_rbegin->time != timeline_begin->time) {
          return absl::InternalError(
              absl::StrCat("Duplicated time mismatch ",
                           DebugString(timeline_rbegin->time),
                           " and ",
                           DebugString(timeline_begin->time),
                           " for segment #",
                           i));
        } else if (timeline_rbegin->degrees_of_freedom !=
                   timeline_begin->degrees_of_freedom) {
          bool const left_nan = timeline_rbegin->degrees_of_freedom !=
                                timeline_rbegin->degrees_of_freedom;
          bool const right_nan = timeline_begin->degrees_of_freedom !=
                                 timeline_begin->degrees_of_freedom;
          if (!(left_nan && right_nan)) {
            return absl::InternalError(
                absl::StrCat("Duplicated degrees of freedom mismatch ",
                             DebugString(timeline_rbegin->degrees_of_freedom),
                             " and ",
                             DebugString(timeline_begin->degrees_of_freedom),
                             " for segment #",
                             i));
          }
        }
      }
    }
  }
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::AdjustAfterSplicing(
    DiscreteTrajectory& from,
    DiscreteTrajectory& to,
    typename Segments::iterator to_segments_begin) {
  // Reset the self pointers of the new segments.  This is necessary for
  // iterating over them.
  for (auto sit = to_segments_begin; sit != to.segments_->end(); ++sit) {
    sit->SetSelf(SegmentIterator(to.segments_.get(), sit));
  }

  // Copy the time-to-segment entries on or after the time of |to_segment_begin|
  // from |from| to |to|.  We are sure to copy all the entries we need because
  // it includes the last segment that starts at that time.
  auto from_leit = from.FindSegment(to_segments_begin->front().time);
  for (auto leit = from_leit;
       leit != from.segment_by_left_endpoint_.end();
       ++leit) {
    to.segment_by_left_endpoint_.insert_or_assign(
        to.segment_by_left_endpoint_.end(), leit->first, leit->second);
  }

  // Erase from |from| the entries that we just copied.  We may erase too much
  // if |from| now ends with a 1-point segment that was not in the map, so we
  // insert it again.
  from.segment_by_left_endpoint_.erase(from_leit,
                                       from.segment_by_left_endpoint_.end());
  if (!from.segments_->empty()) {
    auto const last_segment = --from.segments_->end();
    if (!last_segment->empty()) {
      from.segment_by_left_endpoint_.insert_or_assign(
          from.segment_by_left_endpoint_.end(),
          last_segment->front().time,
          last_segment);
    }
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ReadFromPreHamiltonMessage(
    serialization::DiscreteTrajectory::Downsampling const& message,
    DownsamplingParameters& downsampling_parameters,
    Instant& start_of_dense_timeline) {
  bool const is_pre_haar = message.has_start_of_dense_timeline();
  LOG_IF(WARNING, is_pre_haar)
      << "Reading pre-Haar DiscreteTrajectory.Downsampling";

  downsampling_parameters = {
      .max_dense_intervals = message.max_dense_intervals(),
      .tolerance = Length::ReadFromMessage(message.tolerance())};
  if (is_pre_haar) {
    start_of_dense_timeline =
        Instant::ReadFromMessage(message.start_of_dense_timeline());
  } else {
    start_of_dense_timeline =
        Instant::ReadFromMessage(message.dense_timeline(0));
  }
}

template<typename Frame>
DiscreteTrajectorySegmentIterator<Frame>
DiscreteTrajectory<Frame>::ReadFromPreHamiltonMessage(
    serialization::DiscreteTrajectory::Brood const& message,
    std::vector<SegmentIterator*> const& tracked,
    value_type const& fork_point,
    DiscreteTrajectory& trajectory) {
  CHECK_EQ(fork_point.time, Instant::ReadFromMessage(message.fork_time()))
      << "Cannot read trajectory with a fork not at end of segment";
  CHECK_EQ(1, message.trajectories_size())
      << "Cannot read trajectory with multiple forks";

  // Keep an iterator to the last segment to be able to return a segment
  // iterator.
  auto sit = --trajectory.segments_->end();
  ReadFromPreHamiltonMessage(
      message.trajectories(0), tracked, fork_point, trajectory);
  ++sit;

  return SegmentIterator(trajectory.segments_.get(), sit);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ReadFromPreHamiltonMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked,
    std::optional<value_type> const& fork_point,
    DiscreteTrajectory& trajectory) {
  bool const is_pre_frobenius = !message.has_zfp();
  LOG_IF(WARNING, is_pre_frobenius)
      << "Reading pre-Frobenius DiscreteTrajectory";

  trajectory.segments_->emplace_back();
  auto const sit = --trajectory.segments_->end();
  auto const self = SegmentIterator(trajectory.segments_.get(), sit);
  if (is_pre_frobenius) {
    // Pre-Frobenius saves don't use ZFP so we reconstruct them by appending
    // points to the segment.  Note that this must happens before restoring the
    // downsampling parameters to avoid re-downsampling.
    *sit = DiscreteTrajectorySegment<Frame>(self);
    for (auto const& instantaneous_dof : message.timeline()) {
      sit->Append(Instant::ReadFromMessage(instantaneous_dof.instant()),
                  DegreesOfFreedom<Frame>::ReadFromMessage(
                      instantaneous_dof.degrees_of_freedom())).IgnoreError();
    }
    if (message.has_downsampling()) {
      DownsamplingParameters downsampling_parameters;
      Instant start_of_dense_timeline;
      ReadFromPreHamiltonMessage(message.downsampling(),
                                 downsampling_parameters,
                                 start_of_dense_timeline);
      sit->SetDownsamplingUnconditionally(downsampling_parameters);
      sit->SetStartOfDenseTimeline(start_of_dense_timeline);
    }
  } else {
    // Starting with Frobenius we use ZFP so the easiest is to build a
    // serialized segment from the fields of the legacy message and read from
    // that serialized segment.  Note that restoring the number of dense points
    // can only happen once we have reconstructed the timeline.
    serialization::DiscreteTrajectorySegment serialized_segment;
    *serialized_segment.mutable_zfp() = message.zfp();
    *serialized_segment.mutable_exact() = message.exact();

    DownsamplingParameters downsampling_parameters;
    Instant start_of_dense_timeline;
    if (message.has_downsampling()) {
      ReadFromPreHamiltonMessage(message.downsampling(),
                                 downsampling_parameters,
                                 start_of_dense_timeline);
      auto* const serialized_downsampling_parameters =
          serialized_segment.mutable_downsampling_parameters();
      serialized_downsampling_parameters->set_max_dense_intervals(
          downsampling_parameters.max_dense_intervals);
      downsampling_parameters.tolerance.WriteToMessage(
          serialized_downsampling_parameters->mutable_tolerance());
      serialized_segment.set_number_of_dense_points(0);  // Overridden later.
    }
    *sit = DiscreteTrajectorySegment<Frame>::ReadFromMessage(serialized_segment,
                                                             self);
    if (message.has_downsampling()) {
      sit->SetStartOfDenseTimeline(start_of_dense_timeline);
    }
  }

  // Create the duplicated point if needed.
  if (fork_point.has_value()) {
    sit->SetForkPoint(fork_point.value());
  }

  // Restore the (single) child as the next segment.
  if (message.children_size() == 1) {
    auto const child = ReadFromPreHamiltonMessage(message.children(0),
                                                  tracked,
                                                  /*fork_point=*/*sit->rbegin(),
                                                  trajectory);

    // There were no fork positions prior to Буняковский.
    bool const has_fork_position = message.fork_position_size() > 0;
    if (has_fork_position) {
      CHECK_EQ(1, message.fork_position_size())
          << "Cannot read trajectory with " << message.fork_position_size()
          << " fork positions";
      int const fork_position = message.fork_position(0);
      *tracked[fork_position] = child;
    }
  } else if (message.children_size() > 1) {
    LOG(FATAL) << "Cannot read trajectory with " << message.children_size()
               << " children";
  }

  // Finally, set the time-to-segment map.
  if (!sit->empty()) {
    // This is the *only* place where we must use |emplace|, not
    // |insert_or_assign|.  The reason is that this happens when returning from
    // the recursivity (see to the call to ReadFromPreHamiltonMessage) so
    // segments are processed in reverse order.  Therefore, a segment that is
    // the last at its time will be processed *before* any 1-point segments with
    // the same time, and must be the one stored in the map.
    trajectory.segment_by_left_endpoint_.emplace(sit->begin()->time, sit);
  }
}

}  // namespace internal_discrete_trajectory
}  // namespace physics
}  // namespace principia
