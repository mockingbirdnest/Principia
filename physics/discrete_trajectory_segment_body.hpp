#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include <algorithm>
#include <iterator>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/btree_set.h"
#include "base/zfp_compressor.hpp"
#include "glog/logging.h"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _discrete_trajectory_segment {
namespace internal {

using namespace principia::base::_zfp_compressor;
using namespace principia::quantities::_si;

template<typename Frame>
DiscreteTrajectorySegment<Frame>::DiscreteTrajectorySegment(
    DiscreteTrajectorySegmentIterator<Frame> const self)
    : self_(self) {}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetDownsamplingUnconditionally(
    DownsamplingParameters const& downsampling_parameters) {
  downsampling_parameters_ = downsampling_parameters;
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reference
DiscreteTrajectorySegment<Frame>::front() const {
  return *begin();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reference
DiscreteTrajectorySegment<Frame>::back() const {
  return *std::prev(timeline_.end());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::begin() const {
  return iterator(self_, timeline_.begin());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::end() const {
  if (timeline_.empty()) {
    return iterator(self_, timeline_.end());
  } else {
    // The decrement/increment ensures that we normalize the end iterator to the
    // next segment or to the end of the trajectory.  This is relatively
    // expensive, taking 25-30 ns.
    return ++iterator(self_, --timeline_.end());
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rbegin() const {
  return reverse_iterator(end());
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::reverse_iterator
DiscreteTrajectorySegment<Frame>::rend() const {
  return reverse_iterator(begin());
}

template<typename Frame>
bool DiscreteTrajectorySegment<Frame>::empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::size() const {
  // NOTE(phl): This assumes that there are no repeated times *within* a
  // segment.  This is enforced by Append.
  return timeline_.size();
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::clear() {
  downsampling_parameters_.reset();
  downsampling_error_ = Length{};
  timeline_.clear();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::find(Instant const& t) const {
  auto const it = FindOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::lower_bound(Instant const& t) const {
  auto const it = LowerBoundOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::upper_bound(Instant const& t) const {
  auto const it = UpperBoundOrNullopt(t);
  return it.has_value() ? *it : end();
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_min() const {
  return timeline_.empty() ? InfiniteFuture : timeline_.cbegin()->time;
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_max() const {
  return timeline_.empty() ? InfinitePast : timeline_.crbegin()->time;
}

template<typename Frame>
Position<Frame> DiscreteTrajectorySegment<Frame>::EvaluatePosition(
    Instant const& t) const {
  // Doing the checks first is slightly more efficient than after the lookup.
  CHECK_LE(t_min(), t);
  CHECK_GE(t_max(), t);
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom.position();
  }
  return get_hermite3(it)(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectorySegment<Frame>::EvaluateVelocity(
    Instant const& t) const {
  // Doing the checks first is slightly more efficient than after the lookup.
  CHECK_LE(t_min(), t);
  CHECK_GE(t_max(), t);
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom.velocity();
  }
  return get_hermite3(it).EvaluateDerivative(t);
}

template<typename Frame>
DegreesOfFreedom<Frame>
DiscreteTrajectorySegment<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  // Doing the checks first is slightly more efficient than after the lookup.
  CHECK_LE(t_min(), t);
  CHECK_GE(t_max(), t);
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom;
  }
  auto const& interpolation = get_hermite3(it);
  return {interpolation(t), interpolation.EvaluateDerivative(t)};
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetDownsampling(
    DownsamplingParameters const& downsampling_parameters) {
  // The semantics of changing downsampling on a segment that has 2 points or
  // more are unclear.  Let's not do that.
  CHECK_LE(timeline_.size(), 1);
  downsampling_parameters_ = downsampling_parameters;
  downsampling_error_ = Length{};
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ClearDownsampling() {
  downsampling_parameters_ = std::nullopt;
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectorySegment*> message,
    std::vector<iterator> const& exact) const {
  WriteToMessage(message,
                 timeline_.begin(),
                 timeline_.end(),
                 timeline_.size(),
                 /*number_of_points_to_skip_at_end*/ 0,
                 exact);
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectorySegment*> message,
    iterator const begin,
    iterator const end,
    std::vector<iterator> const& exact) const {
  typename Timeline::const_iterator const timeline_begin =
      begin == this->end() ? timeline_.cend()
                           : iterator::iterator(begin.point_);
  typename Timeline::const_iterator const timeline_end =
      end == this->end() ? timeline_.cend()
                         : iterator::iterator(end.point_);
  bool const covers_entire_segment =
      timeline_begin == timeline_.cbegin() && timeline_end == timeline_.cend();
  std::int64_t const timeline_size =
      covers_entire_segment ? timeline_.size()
                            : std::distance(timeline_begin, timeline_end);
  std::int64_t const number_of_points_to_skip_at_end =
      covers_entire_segment ? 0 : std::distance(timeline_end, timeline_.end());
  WriteToMessage(message,
                 timeline_begin,
                 timeline_end,
                 timeline_size,
                 number_of_points_to_skip_at_end,
                 exact);
}

template<typename Frame>
DiscreteTrajectorySegment<Frame>
DiscreteTrajectorySegment<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectorySegment const& message,
    DiscreteTrajectorySegmentIterator<Frame> const self)
  requires serializable<Frame> {
  bool const is_pre_лефшец = !message.zfp().has_is_lefschetz_timeline();
  LOG_IF(WARNING, is_pre_лефшец)
      << "Reading pre-Лефшец DiscreteTrajectorySegment";

  DiscreteTrajectorySegment<Frame> segment(self);

  // Construct a map for efficient lookup of the exact points.
  Timeline exact;
  for (auto const& instantaneous_degrees_of_freedom : message.exact()) {
    exact.emplace_hint(
        exact.cend(),
        Instant::ReadFromMessage(instantaneous_degrees_of_freedom.instant()),
        DegreesOfFreedom<Frame>::ReadFromMessage(
            instantaneous_degrees_of_freedom.degrees_of_freedom()));
  }

  // Decompress the timeline before restoring the downsampling parameters to
  // avoid re-downsampling.
  ZfpCompressor decompressor;
  ZfpCompressor::ReadVersion(message);

  int const timeline_size = message.zfp().timeline_size();
  std::vector<double> t(timeline_size);
  std::vector<double> qx(timeline_size);
  std::vector<double> qy(timeline_size);
  std::vector<double> qz(timeline_size);
  std::vector<double> px(timeline_size);
  std::vector<double> py(timeline_size);
  std::vector<double> pz(timeline_size);
  std::vector<double> err(timeline_size);
  std::string_view zfp_timeline(message.zfp().timeline().data(),
                                message.zfp().timeline().size());

  decompressor.ReadFromMessageMultidimensional<2>(t, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qx, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qy, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qz, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(px, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(py, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(pz, zfp_timeline);
  if (is_pre_лефшец) {
    // Pre-Лефшец saves didn't store the downsampling error; we have to assume
    // the worst for the points that were downsampled, and no error for the
    // dense points.
    std::int64_t const downsampled_size =
        timeline_size - message.number_of_dense_points();
    for (std::int64_t i = 0; i < downsampled_size; ++i) {
      err[i] = Infinity<double>;
    }
    for (std::int64_t i = downsampled_size; i < timeline_size; ++i) {
      err[i] = 0;
    }
  } else {
    decompressor.ReadFromMessageMultidimensional<2>(err, zfp_timeline);
  }

  for (int i = 0; i < timeline_size; ++i) {
    Position<Frame> const q =
        Frame::origin +
        Displacement<Frame>({qx[i] * Metre, qy[i] * Metre, qz[i] * Metre});
    Velocity<Frame> const p({px[i] * (Metre / Second),
                             py[i] * (Metre / Second),
                             pz[i] * (Metre / Second)});
    Length const downsampling_error = err[i] * Metre;

    // See if this is a point whose degrees of freedom must be restored
    // exactly.  Note that the calls to `Append` recompute the interpolation as
    // being exact since we didn't set the downsampling parameters yet.
    Instant const time = Instant() + t[i] * Second;
    if (auto it = exact.find(time); it == exact.cend()) {
      segment.Append(time, DegreesOfFreedom<Frame>(q, p)).IgnoreError();
    } else {
      segment.Append(time, it->degrees_of_freedom).IgnoreError();
    }

    // Restore the downsampling error for the interpolation ending at this
    // point.  Remember that there is no interpolation for the first point.
    if (i > 0) {
      segment.timeline_.rbegin()->interpolation->error = downsampling_error;
    }
  }

  // Finally, restore the downsampling information.
  if (message.has_downsampling_parameters()) {
    segment.downsampling_parameters_ = DownsamplingParameters{
        .tolerance = Length::ReadFromMessage(
            message.downsampling_parameters().tolerance())};
  }

  return segment;
}

template <typename Frame>
std::optional<typename DiscreteTrajectorySegment<Frame>::iterator>
DiscreteTrajectorySegment<Frame>::FindOrNullopt(Instant const& t) const {
  auto const it = timeline_.find(t);
  if (it == timeline_.end()) {
    return std::nullopt;
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
std::optional<typename DiscreteTrajectorySegment<Frame>::iterator>
DiscreteTrajectorySegment<Frame>::LowerBoundOrNullopt(Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it == timeline_.end()) {
    return std::nullopt;
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
std::optional<typename DiscreteTrajectorySegment<Frame>::iterator>
DiscreteTrajectorySegment<Frame>::UpperBoundOrNullopt(Instant const& t) const {
  auto const it = timeline_.upper_bound(t);
  if (it == timeline_.end()) {
    return std::nullopt;
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetSelf(
    DiscreteTrajectorySegmentIterator<Frame> const self) {
  self_ = self;
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::Prepend(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (!timeline_.empty()) {
    CHECK(t < timeline_.cbegin()->time)
        << "Prepend out of order at " << t << ", first time is "
        << timeline_.cbegin()->time;
  }
  auto const it = timeline_.emplace_hint(timeline_.cbegin(),
                                         t,
                                         degrees_of_freedom);

  if (auto const next = std::next(it); next != timeline_.cend()) {
    CreateInterpolation(next);
  }
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(Instant const& t) {
  ForgetAfter(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(
    typename Timeline::const_iterator const begin) {
  timeline_.erase(begin, timeline_.cend());
  if (!timeline_.empty()) {
    // Restore the downsampling error resulting from the last interpolation, if
    // there is one.
    auto const& last_interpolation = timeline_.crbegin()->interpolation;
    if (last_interpolation != nullptr) {
      downsampling_error_ = last_interpolation->error;
      return;
    }
  }
  downsampling_error_ = Length{};
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(Instant const& t) {
  ForgetBefore(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(
    typename Timeline::const_iterator const end) {
  std::int64_t const number_of_points_to_remove =
      std::distance(timeline_.cbegin(), end);
  timeline_.erase(timeline_.cbegin(), end);
  if (!timeline_.empty()) {
    timeline_.begin()->interpolation = nullptr;
  }
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (timeline_.empty()) {
    // The first point of the segment, no interpolation yet (but this will be
    // the start point for future interpolations).
    auto const it =
        timeline_.emplace_hint(timeline_.cend(), t, degrees_of_freedom);
  } else {
    auto const first = timeline_.cbegin();
    auto const last = std::prev(timeline_.cend());
    if (first->time == t) {
      LOG(WARNING) << "Append at existing time " << t << ", time range = ["
                   << first->time << ", " << last->time << "]";
      return absl::OkStatus();
    }

    CHECK_LT(last->time, t)
        << "Append out of order at " << t << ", last time is "
        << last->time;

    if (timeline_.size() > 1 && downsampling_parameters_.has_value()) {
      // The lower point of the interpolation is the penultimate point of the
      // segment, since we don't retain the intermediate points that were below
      // the tolerance.  Build an interpolation from that point to the new point
      // being added.
      auto const lower = std::prev(last);
      auto hermite3 = MakeHermite3(lower, t, degrees_of_freedom);

      // Compute the new error bound.
      auto const hermite3_difference =
          hermite3 - last->interpolation->hermite3;
      downsampling_error_ += hermite3_difference.LInfinityL₂Norm();
      if (downsampling_error_ < downsampling_parameters_->tolerance) {
        // The error bound is below the tolerance, replace the last point with
        // the new one, record the new interpolation, and keep going.  This is
        // faster than `erase` followed by `emplace_hint`.
        auto extracted = timeline_.extract(last);
        auto& value = extracted.value();
        value.time = t;
        value.degrees_of_freedom = degrees_of_freedom;
        value.interpolation =
            NewInterpolation(std::move(hermite3), downsampling_error_);
        timeline_.insert(timeline_.cend(), std::move(extracted));
        return absl::OkStatus();
      }
    }

    // We land here if (1) the point being appended is the second point of the
    // segment; or (2) the error bound is above the tolerance; or (3)
    // downsampling is disabled.  In all cases, we need to build an
    // interpolation starting at the last point and append
    // our new point.
    auto const lower = last;
    auto hermite3 = MakeHermite3(lower, t, degrees_of_freedom);
    auto const it =
        timeline_.emplace_hint(timeline_.cend(), t, degrees_of_freedom);
    // There is no error because the points being interpolated are consecutive,
    // so the interpolation exactly matches their positions and velocities.
    downsampling_error_ = Length{};
    it->interpolation =
        NewInterpolation(std::move(hermite3), downsampling_error_);
  }

  return absl::OkStatus();
}

// Ideally, the segment constructed by reanimation should end with exactly the
// same time and degrees of freedom as the start of the non-collapsible segment.
// Unfortunately, we believe that numerical inaccuracies are introduced by the
// computations that go through parts, and this introduces small errors.
// TODO(egg): Change Vessel to use PileUp directly and not go through Part.
#define PRINCIPIA_MERGE_STRICT_CONSISTENCY 0

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::Merge(
    DiscreteTrajectorySegment<Frame> segment) {
  if (segment.timeline_.empty()) {
    return;
  } else if (timeline_.empty()) {
    downsampling_parameters_ = segment.downsampling_parameters_;
    timeline_ = std::move(segment.timeline_);
  } else if (auto const [this_crbegin, segment_begin] =
                 std::pair{std::prev(timeline_.cend()),
                           segment.timeline_.begin()};
             this_crbegin->time <= segment_begin->time) {
#if PRINCIPIA_MERGE_STRICT_CONSISTENCY
    CHECK(this_crbegin->time < segment_cbegin->time ||
          this_crbegin->degrees_of_freedom ==
              segment_cbegin->degrees_of_freedom)
        << "Inconsistent merge: [" << segment.timeline_.cbegin()->time
        << ", " << std::prev(segment.timeline_.cend())->time
        << "] into [" << timeline_.cbegin()->time
        << ", " << std::prev(timeline_.cend())->time
        << "], degrees_of_freedom "
        << this_crbegin->degrees_of_freedom << " and "
        << segment_cbegin->degrees_of_freedom << " don't match";
#endif
    Instant const segment_begin_time = segment_begin->time;
    downsampling_parameters_ = segment.downsampling_parameters_;
    timeline_.merge(segment.timeline_);
    // There may not be an interpolation at `segment_begin_time` (there may be
    // one if the segments have a common time, though). (Re)compute it, but
    // remember that we cannot trust that `segment_begin` is in `timeline_` (or
    // even valid) so we need a lookup.
    auto const it = timeline_.find(segment_begin_time);
    CHECK(it != timeline_.cend());
    if (it != timeline_.cbegin()) {
      UpdateInterpolation(it);
    }
  } else if (auto const [segment_crbegin, this_begin] =
                 std::pair{std::prev(segment.timeline_.cend()),
                           timeline_.begin()};
             segment_crbegin->time <= this_begin->time) {
#if PRINCIPIA_MERGE_STRICT_CONSISTENCY
    CHECK(segment_crbegin->time < this_cbegin->time ||
          segment_crbegin->degrees_of_freedom ==
              this_cbegin->degrees_of_freedom)
        << "Inconsistent merge: [" << segment.timeline_.cbegin()->time
        << ", " << std::prev(segment.timeline_.cend())->time
        << "] into [" << timeline_.cbegin()->time
        << ", " << std::prev(timeline_.cend())->time
        << "], degrees_of_freedom "
        << segment_crbegin->degrees_of_freedom << " and "
        << this_cbegin->degrees_of_freedom << " don't match";
#endif
    Instant const this_begin_time = this_begin->time;
    timeline_.merge(segment.timeline_);
    // There may not be an interpolation at `this_begin_time` (there may be one
    // if the segments have a common time, though). (Re)compute it, but remember
    // that we cannot trust that `this_begin` is in `timeline_` (or even valid)
    // so we need a lookup.
    auto const it = timeline_.find(this_begin_time);
    CHECK(it != timeline_.cend());
    if (it != timeline_.cbegin()) {
      UpdateInterpolation(it);
    }
  } else {
    LOG(FATAL) << "Overlapping merge: [" << segment.timeline_.cbegin()->time
               << ", " << std::prev(segment.timeline_.cend())->time
               << "] into [" << timeline_.cbegin()->time
               << ", " << std::prev(timeline_.cend())->time << "]";
  }
}

#undef PRINCIPIA_MERGE_STRICT_CONSISTENCY

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetForkPoint(value_type const& point) {
  auto const it = timeline_.emplace_hint(
      timeline_.begin(), point.time, point.degrees_of_freedom);
  CHECK(it == timeline_.begin())
      << "Inconsistent fork point at time " << point.time;

  if (auto const next = std::next(it); next != timeline_.cend()) {
    CreateInterpolation(next);
  }
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::CreateInterpolation(
    typename Timeline::iterator const upper) {
  CHECK(upper != timeline_.cbegin());
  CHECK_EQ(nullptr, upper->interpolation)
      << "Already has an interpolation at " << upper->time;
  UpdateInterpolation(upper);
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::UpdateInterpolation(
    typename Timeline::iterator const upper) {
  CHECK(upper != timeline_.cbegin());
  auto const lower = std::prev(upper);
  auto const& [lower_time, lower_degrees_of_freedom] = *lower;
  auto const& [upper_time, upper_degrees_of_freedom] = *upper;
  upper->interpolation =
      NewInterpolation(Hermite3<Position<Frame>, Instant>(
                           std::pair{lower_time, upper_time},
                           std::pair{lower_degrees_of_freedom.position(),
                                     upper_degrees_of_freedom.position()},
                           std::pair{lower_degrees_of_freedom.velocity(),
                                     upper_degrees_of_freedom.velocity()}),
                       /*error=*/Length{});
}

template<typename Frame>
Hermite3<Position<Frame>, Instant>
DiscreteTrajectorySegment<Frame>::MakeHermite3(
    typename Timeline::const_iterator lower,
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  auto const lower_time = lower->time;
  auto const lower_degrees_of_freedom = lower->degrees_of_freedom;
  return Hermite3<Position<Frame>, Instant>(
      std::pair{lower_time, t},
      std::pair{lower_degrees_of_freedom.position(),
                degrees_of_freedom.position()},
      std::pair{lower_degrees_of_freedom.velocity(),
                degrees_of_freedom.velocity()});
}

template<typename Frame>
not_null<std::unique_ptr<Interpolation<Frame>>>
DiscreteTrajectorySegment<Frame>::NewInterpolation(
    Hermite3<Position<Frame>, Instant> hermite3,
    Length const& error) {
  return make_not_null_unique<Interpolation<Frame>>(std::move(hermite3), error);
}

template<typename Frame>
Hermite3<Position<Frame>, Instant> const&
DiscreteTrajectorySegment<Frame>::get_hermite3(
    typename Timeline::const_iterator const upper) const {
  CHECK(upper != timeline_.cbegin());
  CHECK_NOTNULL(upper->interpolation);
  return upper->interpolation->hermite3;
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_begin() const {
  return timeline_.cbegin();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::Timeline::const_iterator
DiscreteTrajectorySegment<Frame>::timeline_end() const {
  return timeline_.cend();
}

template<typename Frame>
bool DiscreteTrajectorySegment<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectorySegment<Frame>::timeline_size() const {
  return timeline_.size();
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectorySegment*> message,
    typename Timeline::const_iterator const timeline_begin,
    typename Timeline::const_iterator const timeline_end,
    std::int64_t const timeline_size,
    std::int64_t const number_of_points_to_skip_at_end,
    std::vector<iterator> const& exact) const {
  if (downsampling_parameters_.has_value()) {
    auto* const serialized_downsampling_parameters =
        message->mutable_downsampling_parameters();
    downsampling_parameters_->tolerance.WriteToMessage(
        serialized_downsampling_parameters->mutable_tolerance());
  }

  // Convert the `exact` vector into a set, and add the extremities.  This
  // ensures that we don't have redundancies.  The set is sorted by time to
  // guarantee that serialization is reproducible.
  auto time_comparator = [](value_type const* const left,
                            value_type const* const right) {
    return left->time < right->time;
  };
  absl::btree_set<value_type const*,
                  decltype(time_comparator)> exact_set(time_comparator);
  for (auto const it : exact) {
    exact_set.insert(&*it);
  }
  if (timeline_size > 0) {
    exact_set.insert(&*timeline_begin);
    exact_set.insert(&*std::prev(timeline_end));
  }

  // Serialize the exact points.
  for (auto const* ptr : exact_set) {
    auto const& [t, degrees_of_freedom] = *ptr;
    auto* const serialized_exact = message->add_exact();
    t.WriteToMessage(serialized_exact->mutable_instant());
    degrees_of_freedom.WriteToMessage(
        serialized_exact->mutable_degrees_of_freedom());
  }

  auto* const zfp = message->mutable_zfp();
  zfp->set_timeline_size(timeline_size);
  zfp->set_is_lefschetz_timeline(true);

  // The timeline data is made dimensionless and stored in separate arrays per
  // coordinate.  We expect strong correlations within a coordinate over time,
  // but not between coordinates.
  std::vector<double> t;
  std::vector<double> qx;
  std::vector<double> qy;
  std::vector<double> qz;
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> err;
  t.reserve(timeline_size);
  qx.reserve(timeline_size);
  qy.reserve(timeline_size);
  qz.reserve(timeline_size);
  px.reserve(timeline_size);
  py.reserve(timeline_size);
  pz.reserve(timeline_size);
  err.reserve(timeline_size);
  std::optional<Instant> previous_instant;
  Time max_Δt;
  std::string* const zfp_timeline = zfp->mutable_timeline();
  for (auto it = timeline_begin; it != timeline_end; ++it) {
    auto const& [instant, degrees_of_freedom] = *it;
    auto const q = degrees_of_freedom.position() - Frame::origin;
    auto const p = degrees_of_freedom.velocity();
    t.push_back((instant - Instant{}) / Second);
    qx.push_back(q.coordinates().x / Metre);
    qy.push_back(q.coordinates().y / Metre);
    qz.push_back(q.coordinates().z / Metre);
    px.push_back(p.coordinates().x / (Metre / Second));
    py.push_back(p.coordinates().y / (Metre / Second));
    pz.push_back(p.coordinates().z / (Metre / Second));
    err.push_back(it->interpolation != nullptr
                      ? it->interpolation->error / Metre : 0.0);
    if (previous_instant.has_value()) {
      max_Δt = std::max(max_Δt, instant - *previous_instant);
    }
    previous_instant = instant;
  }

  // Times and errors are exact.
  ZfpCompressor time_compressor(0);
  ZfpCompressor error_compressor(0);
  // Lengths are approximated to the downsampling tolerance if downsampling is
  // enabled, otherwise they are exact.
  Length const length_tolerance = downsampling_parameters_.has_value()
                                      ? downsampling_parameters_->tolerance
                                      : Length();
  ZfpCompressor length_compressor(length_tolerance / Metre);
  // Speeds are approximated based on the length tolerance and the maximum
  // step in the timeline.
  ZfpCompressor const speed_compressor((length_tolerance / max_Δt) /
                                        (Metre / Second));

  ZfpCompressor::WriteVersion(message);
  time_compressor.WriteToMessageMultidimensional<2>(t, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qx, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qy, zfp_timeline);
  length_compressor.WriteToMessageMultidimensional<2>(qz, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(px, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(py, zfp_timeline);
  speed_compressor.WriteToMessageMultidimensional<2>(pz, zfp_timeline);
  error_compressor.WriteToMessageMultidimensional<2>(err, zfp_timeline);
}

}  // namespace internal
}  // namespace _discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
