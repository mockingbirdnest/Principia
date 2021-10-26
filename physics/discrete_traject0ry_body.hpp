#pragma once

#include "physics/discrete_traject0ry.hpp"

#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "base/status_utilities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_traject0ry {

using astronomy::InfinitePast;
using base::make_not_null_unique;
using base::uninitialized;
using quantities::Length;

template<typename Frame>
DiscreteTraject0ry<Frame>::DiscreteTraject0ry()
    : segments_(make_not_null_unique<Segments>(1)) {
  auto const sit = segments_->begin();
  auto const self = SegmentIterator(segments_.get(), sit);
  *sit = DiscreteTrajectorySegment<Frame>(self);
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::reference
DiscreteTraject0ry<Frame>::front() const {
  return *begin();
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::reference
DiscreteTraject0ry<Frame>::back() const {
  return *rbegin();
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::iterator
DiscreteTraject0ry<Frame>::begin() const {
  return segments_->front().begin();
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::iterator
DiscreteTraject0ry<Frame>::end() const {
  return segments_->back().end();
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::reverse_iterator
DiscreteTraject0ry<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::reverse_iterator
DiscreteTraject0ry<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTraject0ry<Frame>::empty() const {
  return segment_by_left_endpoint_.empty();
}

template<typename Frame>
std::int64_t DiscreteTraject0ry<Frame>::size() const {
  std::int64_t size = 1;
  for (auto const& segment : *segments_) {
    size += segment.size();
  }
  size -= segments_->size();  // The junction points.
  return size;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::iterator
DiscreteTraject0ry<Frame>::find(Instant const& t) const {
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
typename DiscreteTraject0ry<Frame>::iterator
DiscreteTraject0ry<Frame>::lower_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    return end();
  }
  auto const sit = leit->second;
  auto const it = sit->lower_bound(t);
  if (it == sit->end()) {
    return end();
  }
  return it;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::iterator
DiscreteTraject0ry<Frame>::upper_bound(Instant const& t) const {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.cend()) {
    return end();
  }
  auto const sit = leit->second;
  auto const it = sit->upper_bound(t);
  if (it == sit->end()) {
    return end();
  }
  return it;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::SegmentRange
DiscreteTraject0ry<Frame>::segments() const {
  return SegmentRange(SegmentIterator(
                          segments_.get(), segments_->begin()),
                      SegmentIterator(
                          segments_.get(), segments_->end()));
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::ReverseSegmentRange
DiscreteTraject0ry<Frame>::rsegments() const {
  return ReverseSegmentRange(std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->end())),
                             std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->begin())));
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::SegmentIterator
DiscreteTraject0ry<Frame>::NewSegment() {
  auto& last_segment = segments_->back();

  auto const& new_segment = segments_->emplace_back();
  auto const new_segment_sit = --segments_->end();
  auto const new_self = SegmentIterator(segments_.get(), new_segment_sit);
  *new_segment_sit = DiscreteTrajectorySegment<Frame>(new_self);

  // It is only possible to insert a segment after an empty segment if the
  // entire trajectory is empty.
  if (last_segment.empty()) {
    CHECK(segment_by_left_endpoint_.empty())
        << "Time to segment map has " << segment_by_left_endpoint_.size()
        << " elements";
  } else {
    auto const& [last_time, last_degrees_of_freedom] = *last_segment.rbegin();
    new_segment_sit->Append(last_time, last_degrees_of_freedom);
    // The use of |insert_or_assign| ensure that we override any entry with the
    // same left endpoint.
    segment_by_left_endpoint_.insert_or_assign(
        segment_by_left_endpoint_.end(), last_time, new_segment_sit);
  }

  return new_self;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::DiscreteTraject0ry
DiscreteTraject0ry<Frame>::DetachSegments(SegmentIterator const begin) {
  DiscreteTraject0ry detached(uninitialized);

  // Move the detached segments to the new trajectory.
  detached.segments_->splice(detached.segments_->end(),
                             *segments_,
                             begin.iterator(), segments_->end());

  AdjustAfterSplicing(/*from=*/*this,
                      /*to=*/detached,
                      /*to_segments_begin=*/detached.segments_->begin());

  return detached;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::SegmentIterator
DiscreteTraject0ry<Frame>::AttachSegments(
    DiscreteTraject0ry&& trajectory) {
  CHECK(!trajectory.empty());
  // NOTE(phl): This check might be too strict, we might want to allow LT as the
  // time comparison, and to adjust the first point of trajectory as needed.
  // We'll see if the clients need that.
  CHECK_EQ(rbegin()->time, trajectory.begin()->time)
      << "Mismatching times when attaching segments";
  CHECK_EQ(rbegin()->degrees_of_freedom, trajectory.begin()->degrees_of_freedom)
      << "Mismatching degrees of freedom when attaching segments";

  if (empty()) {
    *this = DiscreteTraject0ry(uninitialized);
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

  return SegmentIterator(segments_.get(), end_before_splice);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::DeleteSegments(SegmentIterator const begin) {
  segments_->erase(begin.iterator(), segments_->end());
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ForgetAfter(Instant const& t) {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.end()) {
    return;
  }
  auto const sit = leit->second;
  sit->ForgetAfter(t);
  // Here |sit| designates a segment starting at or after |t|.  If |t| is
  // exactly at the beginning of the segment,
  // |DiscreteTrajectorySegment::ForgetAfter| will leave it empty.  In that
  // case we drop the segment entirely.
  if (sit->empty()) {
    segments_->erase(sit, segments_->end());
    segment_by_left_endpoint_.erase(leit, segment_by_left_endpoint_.end());
  } else {
    segments_->erase(std::next(sit), segments_->end());
    segment_by_left_endpoint_.erase(std::next(leit),
                                     segment_by_left_endpoint_.end());
  }
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ForgetAfter(iterator const it) {
  ForgetAfter(it->time);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ForgetBefore(Instant const& t) {
  auto const leit = FindSegment(t);
  if (leit == segment_by_left_endpoint_.end()) {
    return;
  }
  auto const sit = leit->second;
  // This call never makes the segment |*sit| empty.
  sit->ForgetBefore(t);
  for (auto s = segments_->begin(); s != sit; ++s) {
    // This call may either make the segment |*s| empty or leave it with a
    // single point matching |sit->front()|.
    s->ForgetBefore(t);
  }

  // Erase all the entries before and including |leit|.  These are entries for
  // now-empty segments or for 1-point segments, and the segment which we may
  // just have truncated.
  segment_by_left_endpoint_.erase(segment_by_left_endpoint_.begin(),
                                  std::next(leit));
  // Recreate an entry for |sit| with its new left endpoint.
  segment_by_left_endpoint_.insert_or_assign(segment_by_left_endpoint_.begin(),
                                             sit->front().time,
                                             sit);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ForgetBefore(iterator const it) {
  ForgetBefore(it->time);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto leit = FindSegment(t);
  typename Segments::iterator sit;
  if (leit == segment_by_left_endpoint_.end()) {
    // If this is the first point appended to this trajectory, insert it in the
    // time-to-segment map.
    sit = --segments_->end();
    leit = segment_by_left_endpoint_.insert_or_assign(
        segment_by_left_endpoint_.end(), t, sit);
  } else {
    // The segment is expected to always have a point copied from its
    // predecessor.
    sit = leit->second;
    CHECK(!sit->empty()) << "Empty segment at " << t;
  }
  sit->Append(t, degrees_of_freedom);
}

template<typename Frame>
Instant DiscreteTraject0ry<Frame>::t_min() const {
  return segments_->front().t_min();
}

template<typename Frame>
Instant DiscreteTraject0ry<Frame>::t_max() const {
  return segments_->back().t_max();
}

template<typename Frame>
Position<Frame> DiscreteTraject0ry<Frame>::EvaluatePosition(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTraject0ry<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTraject0ry<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t)->second->EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> message,
    std::vector<SegmentIterator> const& tracked,
    std::vector<iterator> const& exact) const {
  // Construct a map to efficiently find if a segment must be tracked.  The
  // keys are pointers to segments in |tracked|, the values are the
  // corresponding indices.
  absl::flat_hash_map<DiscreteTrajectorySegment<Frame> const*, int>
      segment_to_position;
  for (int i = 0; i < tracked.size(); ++i) {
    segment_to_position.emplace(&*tracked[i], i);
  }

  // Initialize the tracked positions to be able to recognize if some are
  // missing.
  message->mutable_tracked_position()->Resize(
      tracked.size(),
      serialization::DiscreteTrajectory::MISSING_TRACKED_POSITION);

  // The position of a segment in the repeated field |segment|.
  int segment_position = 0;
  for (auto sit = segments_->begin();
       sit != segments_->end();
       ++sit, ++segment_position) {
    sit->WriteToMessage(message->add_segment(), exact);

    if (auto const position_it = segment_to_position.find(&*sit);
        position_it != segment_to_position.end()) {
      // The field |tracked_position| is indexed by the indices in |tracked|.
      // Its value is the position of a tracked segment in the field |segment|.
      message->set_tracked_position(position_it->second, segment_position);
    }
  }

  for (auto const& [t, _] : segment_by_left_endpoint_) {
    t.WriteToMessage(message->add_left_endpoint());
  }

  // Check that all the segments in |tracked| were mapped.
  // NOTE(phl): This might be too strong a constraint in Entwurf.
  for (auto const tracked_position : message->tracked_position()) {
    CHECK_NE(serialization::DiscreteTrajectory::MISSING_TRACKED_POSITION,
             tracked_position);
  }
}

template<typename Frame>
template<typename F, typename>
DiscreteTraject0ry<Frame>
DiscreteTraject0ry<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked) {
  DiscreteTraject0ry trajectory(uninitialized);

  bool const is_pre_ζήνων = message.segment_size() == 0;
  if (is_pre_ζήνων) {
    LOG_IF(WARNING, is_pre_ζήνων) << "Reading pre-Ζήνων DiscreteTrajectory";
    ReadFromPreΖήνωνMessage(
        message, tracked, /*fork_point=*/std::nullopt, trajectory);
    CHECK_OK(trajectory.ValidateConsistency());
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
    CHECK_NE(serialization::DiscreteTrajectory::MISSING_TRACKED_POSITION,
             tracked_position);
    *tracked[i] = segment_iterators[tracked_position];
  }

  // Finally restore the left endpoints.
  auto sit = trajectory.segments_->begin();
  for (auto const& serialized_t : message.left_endpoint()) {
    auto const t = Instant::ReadFromMessage(serialized_t);
    trajectory.segment_by_left_endpoint_.insert_or_assign(
        trajectory.segment_by_left_endpoint_.end(), t, sit);
    ++sit;
  }

  CHECK_OK(trajectory.ValidateConsistency());
  return trajectory;
}

template<typename Frame>
DiscreteTraject0ry<Frame>::DiscreteTraject0ry(uninitialized_t)
    : segments_(make_not_null_unique<Segments>()) {}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::SegmentByLeftEndpoint::iterator
DiscreteTraject0ry<Frame>::FindSegment(
    Instant const& t) {
  if (segment_by_left_endpoint_.empty()) {
    return segment_by_left_endpoint_.end();
  }
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin())
      << "No segment covering " << t;
  return --it;
}

template<typename Frame>
typename DiscreteTraject0ry<Frame>::SegmentByLeftEndpoint::const_iterator
DiscreteTraject0ry<Frame>::FindSegment(
    Instant const& t) const {
  if (segment_by_left_endpoint_.empty()) {
    return segment_by_left_endpoint_.cend();
  }
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin())
      << "No segment covering " << t;
  return --it;
}

template<typename Frame>
absl::Status DiscreteTraject0ry<Frame>::ValidateConsistency() const {
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
        absl::StrCat("Times mismatch ",
                      DebugString(left_endpoint),
                      " and ",
                      DebugString(sit->front().time),
                      " between segment #", i, " and the time-to-segment map");
      }
    }
  }
  {
    int i = 0;
    auto sit1 = segments_->cbegin();
    for (auto leit = segment_by_left_endpoint_.cbegin();
         leit != segment_by_left_endpoint_.cend();
         ++leit, ++i) {
      auto const sit2 = leit->second;
      if (sit2->empty()) {
        return absl::InternalError(
            absl::StrCat("Empty segment for time-to-segment entry #", i));
      } else if (sit1 == sit2) {
        ++sit1;
        ++leit;
      } else if (sit1->empty() || sit1->size() == 1) {
        ++sit1;
      } else {
        return absl::InternalError(
            absl::StrCat("Mismatch for time-to-segment entry #", i,
                         " at times ", DebugString(sit1->front().time),
                         " and ", DebugString(sit2->front().time)));
      }
    }
  }
  {
    int i = 0;
    for (auto sit = segments_->cbegin();
         sit != std::prev(segments_->cend());
         ++sit, ++i) {
      // Great care is required here because the DiscreteTrajectoryIterator will
      // "helpfully" paper over differences in degrees of freedom as long as the
      // times match.  We must look at the endpoints of the timeline explicitly.
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
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::AdjustAfterSplicing(
    DiscreteTraject0ry& from,
    DiscreteTraject0ry& to,
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
  if (!from.empty()) {
    auto const last_segment = --from.segments_->end();
    from.segment_by_left_endpoint_.insert_or_assign(
        from.segment_by_left_endpoint_.end(),
        last_segment->front().time,
        last_segment);
  }

}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ReadFromPreΖήνωνMessage(
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
DiscreteTraject0ry<Frame>::ReadFromPreΖήνωνMessage(
    serialization::DiscreteTrajectory::Brood const& message,
    std::vector<SegmentIterator*> const& tracked,
    value_type const& fork_point,
    DiscreteTraject0ry& trajectory) {
  CHECK_EQ(fork_point.time, Instant::ReadFromMessage(message.fork_time()))
      << "Cannot read trajectory with a fork not at end of segment";
  CHECK_EQ(1, message.trajectories_size())
      << "Cannot read trajectory with multiple forks";

  // Keep an iterator to the last segment to be able to return a segment
  // iterator.
  auto sit = --trajectory.segments_->end();
  ReadFromPreΖήνωνMessage(
      message.trajectories(0), tracked, fork_point, trajectory);
  ++sit;

  return SegmentIterator(trajectory.segments_.get(), sit);
}

template<typename Frame>
void DiscreteTraject0ry<Frame>::ReadFromPreΖήνωνMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked,
    std::optional<value_type> const& fork_point,
    DiscreteTraject0ry& trajectory) {
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
                      instantaneous_dof.degrees_of_freedom()));
    }
    if (message.has_downsampling()) {
      DownsamplingParameters downsampling_parameters;
      Instant start_of_dense_timeline;
      ReadFromPreΖήνωνMessage(message.downsampling(),
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
      ReadFromPreΖήνωνMessage(message.downsampling(),
                              downsampling_parameters,
                              start_of_dense_timeline);
      auto* const serialized_downsampling_parameters =
          serialized_segment.mutable_downsampling_parameters();
      serialized_downsampling_parameters->set_max_dense_intervals(
          downsampling_parameters.max_dense_intervals);
      downsampling_parameters.tolerance.WriteToMessage(
          serialized_downsampling_parameters->mutable_tolerance());
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
    auto const child =
        ReadFromPreΖήνωνMessage(message.children(0),
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
  //TODO(phl):Incorrect?
  trajectory.segment_by_left_endpoint_.insert_or_assign(sit->begin()->time,
                                                        sit);
}

}  // namespace internal_discrete_traject0ry
}  // namespace physics
}  // namespace principia
