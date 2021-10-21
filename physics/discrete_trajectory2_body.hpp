#pragma once

#include "physics/discrete_trajectory2.hpp"

#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "base/status_utilities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory2 {

using astronomy::InfinitePast;
using base::make_not_null_unique;
using base::uninitialized;
using quantities::Length;

template<typename Frame>
DiscreteTrajectory2<Frame>::DiscreteTrajectory2()
    : segments_(make_not_null_unique<Segments>(1)) {
  auto const sit = segments_->begin();
  *sit =
      DiscreteTrajectorySegment<Frame>(SegmentIterator(segments_.get(), sit));
  segment_by_left_endpoint_.emplace(InfinitePast, sit);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::begin() const {
  return segments_->front().begin();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::end() const {
  return segments_->back().end();
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::reverse_iterator
DiscreteTrajectory2<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::reverse_iterator
DiscreteTrajectory2<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectory2<Frame>::empty() const {
  for (auto const& segment : *segments_) {
    if (!segment.empty()) {
      return false;
    }
  }
  return true;
}

template<typename Frame>
std::int64_t DiscreteTrajectory2<Frame>::size() const {
  std::int64_t size = 1;
  for (auto const& segment : *segments_) {
    size += segment.size();
  }
  size -= segments_->size();  // The junction points.
  return size;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::find(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->find(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::lower_bound(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->lower_bound(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::iterator
DiscreteTrajectory2<Frame>::upper_bound(Instant const& t) const {
  auto const sit = FindSegment(t);
  auto const it = sit->upper_bound(t);
  if (it == sit->end()) {
    return end();
  } else {
    return it;
  }
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentRange
DiscreteTrajectory2<Frame>::segments() const {
  return SegmentRange(SegmentIterator(
                          segments_.get(), segments_->begin()),
                      SegmentIterator(
                          segments_.get(), segments_->end()));
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::ReverseSegmentRange
DiscreteTrajectory2<Frame>::rsegments() const {
  return ReverseSegmentRange(std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->end())),
                             std::reverse_iterator(SegmentIterator(
                                 segments_.get(), segments_->begin())));
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::NewSegment() {
  auto& last_segment = segments_->back();
  CHECK(!last_segment.empty())
      << "Cannot create a new segment after an empty one";

  auto const& new_segment = segments_->emplace_back();
  auto const new_segment_sit = --segments_->end();
  *new_segment_sit = DiscreteTrajectorySegment<Frame>(
      SegmentIterator(segments_.get(), new_segment_sit));

  auto const& [last_time, last_degrees_of_freedom] = *last_segment.rbegin();
  new_segment_sit->Append(last_time, last_degrees_of_freedom);
  segment_by_left_endpoint_.emplace_hint(
      segment_by_left_endpoint_.end(), last_time, new_segment_sit);

  return SegmentIterator(segments_.get(), new_segment_sit);
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::DiscreteTrajectory2
DiscreteTrajectory2<Frame>::DetachSegments(SegmentIterator const begin) {
  DiscreteTrajectory2 detached(uninitialized);

  // Move the detached segments to the new trajectory.
  detached.segments_->splice(detached.segments_->end(),
                             *segments_,
                             begin.iterator(), segments_->end());

  AdjustAfterSplicing(/*from=*/*this,
                      /*to=*/detached,
                      /*to_segments_begin=*/detached.segments_->begin(),
                      /*to_segments_rend=*/detached.segments_->rend());

  return detached;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::SegmentIterator
DiscreteTrajectory2<Frame>::AttachSegments(
    DiscreteTrajectory2&& trajectory) {
  CHECK(!trajectory.empty());
  // NOTE(phl): This check might be too strict, we might want to allow LT as the
  // time comparison, and to adjust the first point of trajectory as needed.
  // We'll see if the clients need that.
  CHECK_EQ(rbegin()->first, trajectory.begin()->first)
      << "Mismatching times when attaching segments";
  CHECK_EQ(rbegin()->second, trajectory.begin()->second)
      << "Mismatching degrees of freedom when attaching segments";

  if (empty()) {
    *this = DiscreteTrajectory2(uninitialized);
  }

  // The |end| iterator keeps pointing at the end after the splice.  Instead,
  // we track the iterator to the last segment.
  auto const last_before_splice = --segments_->end();

  // Move the attached segments to the end of this trajectory.
  segments_->splice(segments_->end(),
                    *trajectory.segments_);

  auto const end_before_splice = std::next(last_before_splice);
  auto const rbegin_before_splice = std::reverse_iterator(end_before_splice);

  AdjustAfterSplicing(/*from=*/trajectory,
                      /*to=*/*this,
                      /*to_segments_begin=*/end_before_splice,
                      /*to_segments_rend=*/rbegin_before_splice);

  return SegmentIterator(segments_.get(), end_before_splice);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::DeleteSegments(SegmentIterator const begin) {
  segments_->erase(begin.iterator(), segments_->end());
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetAfter(Instant const& t) {
  auto const sit = FindSegment(t);
  sit->ForgetAfter(t);
  // Here |sit| designates a segment starting at or after |t|.  If |t| is
  // exactly at the beginning of the segment, |ForgetAfter| above will leave it
  // empty.  In that case we drop the segment entirely.  Note that this
  // situation doesn't arise for |ForgetBefore| because of the way |FindSegment|
  // works.
  if (sit->empty()) {
    segments_->erase(sit, segments_->end());
  } else {
    segments_->erase(std::next(sit), segments_->end());
  }
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ForgetBefore(Instant const& t) {
  auto const sit = FindSegment(t);
  sit->ForgetBefore(t);
  segments_->erase(segments_->begin(), sit);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto sit = FindSegment(t);
  // If this is the first point appended to this trajectory, insert a proper
  // left endpoint and remove the sentinel.
  if (empty()) {
    segment_by_left_endpoint_.emplace_hint(
        segment_by_left_endpoint_.end(), t, sit);
    segment_by_left_endpoint_.erase(segment_by_left_endpoint_.begin());
  }
  sit->Append(t, degrees_of_freedom);
}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_min() const {
  return segments_->front().t_min();
}

template<typename Frame>
Instant DiscreteTrajectory2<Frame>::t_max() const {
  return segments_->back().t_max();
}

template<typename Frame>
Position<Frame> DiscreteTrajectory2<Frame>::EvaluatePosition(
    Instant const& t) const {
  return FindSegment(t)->EvaluatePosition(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory2<Frame>::EvaluateVelocity(
    Instant const& t) const {
  return FindSegment(t)->EvaluateVelocity(t);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory2<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  return FindSegment(t)->EvaluateDegreesOfFreedom(t);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::WriteToMessage(
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
DiscreteTrajectory2<Frame>
DiscreteTrajectory2<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked) {
  DiscreteTrajectory2 trajectory(uninitialized);

  bool const is_pre_ζήνων = message.segment_size() == 0;
  if (is_pre_ζήνων) {
    LOG_IF(WARNING, is_pre_ζήνων) << "Reading pre-Ζήνων DiscreteTrajectory";
    ReadFromPreΖήνωνMessage(message, tracked, trajectory);
    CHECK_OK(trajectory.IsConsistent());
    return trajectory;
  }

  // First restore the segments themselves.  This vector will be used to restore
  // the tracked segments.
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
    trajectory.segment_by_left_endpoint_.emplace_hint(
        trajectory.segment_by_left_endpoint_.end(), t, sit);
    ++sit;
  }

  CHECK_OK(trajectory.IsConsistent());
  return trajectory;
}

template<typename Frame>
DiscreteTrajectory2<Frame>::DiscreteTrajectory2(uninitialized_t)
    : segments_(make_not_null_unique<Segments>()) {}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::Segments::iterator
DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin()) << "No segment covering " << t;
  return (--it)->second;
}

template<typename Frame>
typename DiscreteTrajectory2<Frame>::Segments::const_iterator
DiscreteTrajectory2<Frame>::FindSegment(
    Instant const& t) const {
  auto it = segment_by_left_endpoint_.upper_bound(t);
  CHECK(it != segment_by_left_endpoint_.begin()) << "No segment covering " << t;
  return (--it)->second;
}

template<typename Frame>
absl::Status DiscreteTrajectory2<Frame>::IsConsistent() const {
  if (segments_->size() != segment_by_left_endpoint_.size()) {
    return absl::InternalError(absl::StrCat("Size mismatch ",
                                            segments_->size(),
                                            " and ",
                                            segment_by_left_endpoint_.size()));
  }
  {
    int i = 0;
    auto it1 = segments_->cbegin();
    auto it2 = segment_by_left_endpoint_.cbegin();
    for (; it1 != segments_->cend(); ++it1, ++it2, ++i) {
      if (it1->begin()->first != it2->first) {
        return absl::InternalError(
            absl::StrCat("Times mismatch ",
                         DebugString(it1->begin()->first),
                         " and ",
                         DebugString(it2->first),
                         " for segment #",
                         i));
      } else if (it1 != it2->second) {
        return absl::InternalError(absl::StrCat("Iterator mismatch ",
                                                " for segment #",
                                                i));
      }
    }
  }
  {
    int i = 0;
    for (auto it = segments_->cbegin();
         it != std::prev(segments_->cend());
         ++it, ++i) {
      if (it->rbegin()->first != std::next(it)->begin()->first) {
        return absl::InternalError(
            absl::StrCat("Duplicated time mismatch ",
                         DebugString(it->rbegin()->first),
                         " and ",
                         DebugString(std::next(it)->begin()->first),
                         " for segment #",
                         i));
      } else if (it->rbegin()->second != std::next(it)->begin()->second) {
        return absl::InternalError(
            absl::StrCat("Duplicated degrees of freedom mismatch ",
                         DebugString(it->rbegin()->second),
                         " and ",
                         DebugString(std::next(it)->begin()->second),
                         " for segment #",
                         i));
      }
    }
  }
  return absl::OkStatus();
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::AdjustAfterSplicing(
    DiscreteTrajectory2& from,
    DiscreteTrajectory2& to,
    typename Segments::iterator to_segments_begin,
    std::reverse_iterator<typename Segments::iterator> to_segments_rend) {

  // Iterate through the target segments to move the time-to-segment mapping
  // from |from| to |to|.
  auto endpoint_it = from.segment_by_left_endpoint_.end();
  for (auto sit = to.segments_->rbegin(); sit != to_segments_rend; ++sit) {
    --endpoint_it;
    to.segment_by_left_endpoint_.insert(
        from.segment_by_left_endpoint_.extract(endpoint_it));
  }

  // Reset the self pointers of the new segments.
  for (auto sit = to_segments_begin; sit != to.segments_->end(); ++sit) {
    sit->SetSelf(SegmentIterator(to.segments_.get(), sit));
  }
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ReadFromPreΖήνωνMessage(
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
DiscreteTrajectory2<Frame>::ReadFromPreΖήνωνMessage(
    serialization::DiscreteTrajectory::Brood const& message,
    std::vector<SegmentIterator*> const& tracked,
    value_type const& fork_point,
    DiscreteTrajectory2& trajectory) {
  CHECK_EQ(fork_point.first, Instant::ReadFromMessage(message.fork_time()))
      << "Cannot read trajectory with a fork not at end of segment";
  CHECK_EQ(1, message.trajectories_size())
      << "Cannot read trajectory with multiple forks";

  // Keep an iterator to the last segment to be able to set the replicated
  // point.
  auto sit = --trajectory.segments_->end();
  ReadFromPreΖήνωνMessage(message.trajectories(0), tracked, trajectory);
  ++sit;
  sit->SetForkPoint(fork_point);

  return SegmentIterator(trajectory.segments_.get(), sit);
}

template<typename Frame>
void DiscreteTrajectory2<Frame>::ReadFromPreΖήνωνMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<SegmentIterator*> const& tracked,
    DiscreteTrajectory2& trajectory) {
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

  // Restore the (single) child as the next segment.
  if (message.children_size() == 1) {
    auto const child = ReadFromPreΖήνωνMessage(message.children(0),
                                               tracked,
                                               *sit->rbegin(),
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
  trajectory.segment_by_left_endpoint_.emplace(sit->begin()->first, sit);
}

}  // namespace internal_discrete_trajectory2
}  // namespace physics
}  // namespace principia
