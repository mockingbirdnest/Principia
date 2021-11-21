#pragma once

#include "physics/discrete_trajectory_segment.hpp"

#include <algorithm>
#include <iterator>
#include <list>
#include <string>
#include <vector>

#include "absl/container/btree_set.h"
#include "astronomy/epoch.hpp"
#include "base/zfp_compressor.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "numerics/fit_hermite_spline.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_discrete_trajectory_segment {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using base::ZfpCompressor;
using geometry::Displacement;
using geometry::Position;
using numerics::FitHermiteSpline;
using quantities::Length;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Second;

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
  return *timeline_.begin();
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
  number_of_dense_points_ = 0;
  timeline_.clear();
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::find(Instant const& t) const {
  auto const it = timeline_.find(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::lower_bound(Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
typename DiscreteTrajectorySegment<Frame>::iterator
DiscreteTrajectorySegment<Frame>::upper_bound(Instant const& t) const {
  auto const it = timeline_.upper_bound(t);
  if (it == timeline_.end()) {
    return end();
  } else {
    return iterator(self_, it);
  }
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_min() const {
  return empty() ? InfiniteFuture : timeline_.cbegin()->time;
}

template<typename Frame>
Instant DiscreteTrajectorySegment<Frame>::t_max() const {
  return empty() ? InfinitePast : timeline_.crbegin()->time;
}

template<typename Frame>
Position<Frame> DiscreteTrajectorySegment<Frame>::EvaluatePosition(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom.position();
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  return GetInterpolation(it).Evaluate(t);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectorySegment<Frame>::EvaluateVelocity(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom.velocity();
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  return GetInterpolation(it).EvaluateDerivative(t);
}

template<typename Frame>
DegreesOfFreedom<Frame>
DiscreteTrajectorySegment<Frame>::EvaluateDegreesOfFreedom(
    Instant const& t) const {
  auto const it = timeline_.lower_bound(t);
  if (it->time == t) {
    return it->degrees_of_freedom;
  }
  CHECK_LT(t_min(), t);
  CHECK_GT(t_max(), t);
  auto const interpolation = GetInterpolation(it);
  return {interpolation.Evaluate(t), interpolation.EvaluateDerivative(t)};
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetDownsampling(
    DownsamplingParameters const& downsampling_parameters) {
  // The semantics of changing downsampling on a segment that has 2 points or
  // more are unclear.  Let's not do that.
  CHECK_LE(timeline_.size(), 1);
  downsampling_parameters_ = downsampling_parameters;
  number_of_dense_points_ = timeline_.empty() ? 0 : 1;
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
template<typename F, typename>
DiscreteTrajectorySegment<Frame>
DiscreteTrajectorySegment<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectorySegment const& message,
    DiscreteTrajectorySegmentIterator<Frame> const self) {
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
  std::string_view zfp_timeline(message.zfp().timeline().data(),
                                message.zfp().timeline().size());

  decompressor.ReadFromMessageMultidimensional<2>(t, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qx, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qy, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(qz, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(px, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(py, zfp_timeline);
  decompressor.ReadFromMessageMultidimensional<2>(pz, zfp_timeline);

  for (int i = 0; i < timeline_size; ++i) {
    Position<Frame> const q =
        Frame::origin +
        Displacement<Frame>({qx[i] * Metre, qy[i] * Metre, qz[i] * Metre});
    Velocity<Frame> const p({px[i] * (Metre / Second),
                             py[i] * (Metre / Second),
                             pz[i] * (Metre / Second)});

    // See if this is a point whose degrees of freedom must be restored
    // exactly.
    Instant const time = Instant() + t[i] * Second;
    if (auto it = exact.find(time); it == exact.cend()) {
      segment.Append(time, DegreesOfFreedom<Frame>(q, p)).IgnoreError();
    } else {
      segment.Append(time, it->degrees_of_freedom).IgnoreError();
    }
  }

  // Finally, restore the downsampling information.
  if (message.has_downsampling_parameters()) {
    segment.downsampling_parameters_ = DownsamplingParameters{
        .max_dense_intervals =
            message.downsampling_parameters().max_dense_intervals(),
        .tolerance = Length::ReadFromMessage(
            message.downsampling_parameters().tolerance())};
  }
  segment.number_of_dense_points_ = message.number_of_dense_points();

  return segment;
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
  CHECK(!timeline_.empty() || t < timeline_.cbegin()->time)
      << "Prepend out of order at " << t << ", first time is "
      << timeline_.cbegin()->time;
  timeline_.emplace_hint(timeline_.cbegin(), t, degrees_of_freedom);
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::Append(
    Instant const& t,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  if (!timeline_.empty() && timeline_.cbegin()->time == t) {
    LOG(WARNING) << "Append at existing time " << t << ", time range = ["
                 << timeline_.cbegin()->time << ", "
                 << timeline_.crbegin()->time << "]";
    return absl::OkStatus();
  }
  auto it = timeline_.emplace_hint(timeline_.cend(),
                                   t,
                                   degrees_of_freedom);
  CHECK(++it == timeline_.end())
      << "Append out of order at " << t << ", last time is "
      << timeline_.crbegin()->time;

  if (downsampling_parameters_.has_value()) {
    return DownsampleIfNeeded();
  } else {
    return absl::OkStatus();
  }
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(Instant const& t) {
  ForgetAfter(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetAfter(
    typename Timeline::const_iterator const begin) {
  std::int64_t number_of_points_to_remove =
      std::distance(begin, timeline_.cend());
  number_of_dense_points_ =
      std::max<std::int64_t>(
          0, number_of_dense_points_ - number_of_points_to_remove);

  timeline_.erase(begin, timeline_.cend());
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(Instant const& t) {
  ForgetBefore(timeline_.lower_bound(t));
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::ForgetBefore(
    typename Timeline::const_iterator const end) {
  std::int64_t number_of_points_to_remove =
      std::distance(timeline_.cbegin(), end);
  number_of_dense_points_ =
      std::max<std::int64_t>(
          0, number_of_dense_points_ - number_of_points_to_remove);

  timeline_.erase(timeline_.cbegin(), end);
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetStartOfDenseTimeline(
    Instant const& t) {
  auto const it = find(t);
  CHECK(it != end()) << "Cannot find time " << t << " in timeline";
  number_of_dense_points_ = std::distance(it, end());
}

template<typename Frame>
void DiscreteTrajectorySegment<Frame>::SetForkPoint(value_type const& point) {
  auto const it = timeline_.emplace_hint(
      timeline_.begin(), point.time, point.degrees_of_freedom);
  CHECK(it == timeline_.begin())
      << "Inconsistent fork point at time " << point.time;
}

template<typename Frame>
absl::Status DiscreteTrajectorySegment<Frame>::DownsampleIfNeeded() {
  ++number_of_dense_points_;
  // Points, hence one more than intervals.
  if (number_of_dense_points_ >
      downsampling_parameters_->max_dense_intervals) {
    // Obtain iterators for all the dense points of the segment.
    using ConstIterators = std::vector<typename Timeline::const_iterator>;
    ConstIterators dense_iterators(number_of_dense_points_);
    CHECK_LE(dense_iterators.size(), timeline_.size());
    auto it = timeline_.crbegin();
    for (int i = dense_iterators.size() - 1; i >= 0; --i) {
      dense_iterators[i] = std::prev(it.base());
      ++it;
    }

    absl::StatusOr<std::list<typename ConstIterators::const_iterator>>
        right_endpoints = FitHermiteSpline<Instant, Position<Frame>>(
            dense_iterators,
            [](auto&& it) -> auto&& { return it->time; },
            [](auto&& it) -> auto&& {
              return it->degrees_of_freedom.position();
            },
            [](auto&& it) -> auto&& {
              return it->degrees_of_freedom.velocity();
            },
            downsampling_parameters_->tolerance);
    if (!right_endpoints.ok()) {
      // Note that the actual appending took place; the propagated status only
      // reflects a lack of downsampling.
      return right_endpoints.status();
    }

    if (right_endpoints->empty()) {
      right_endpoints->push_back(std::prev(dense_iterators.end()));
    }

    // Obtain the times for the right endpoints.  This is necessary because we
    // cannot use iterators for erasing points, as they would get invalidated
    // after the first erasure.
    std::vector<Instant> right_endpoints_times;
    right_endpoints_times.reserve(right_endpoints.value().size());
    for (auto const& it_in_dense_iterators : right_endpoints.value()) {
      right_endpoints_times.push_back((*it_in_dense_iterators)->time);
    }

    // Poke holes in the timeline at the places given by
    // |right_endpoints_times|.  This requires one lookup per erasure.
    auto left_it = dense_iterators.front();
    for (Instant const& right : right_endpoints_times) {
      ++left_it;
      auto const right_it = timeline_.find(right);
      left_it = timeline_.erase(left_it, right_it);
    }
    number_of_dense_points_ = std::distance(left_it, timeline_.cend());
  }
  return absl::OkStatus();
}

template<typename Frame>
Hermite3<Instant, Position<Frame>>
DiscreteTrajectorySegment<Frame>::GetInterpolation(
    typename Timeline::const_iterator const upper) const {
  CHECK(upper != timeline_.cbegin());
  auto const lower = std::prev(upper);
  auto const& [lower_time, lower_degrees_of_freedom] = *lower;
  auto const& [upper_time, upper_degrees_of_freedom] = *upper;
  return Hermite3<Instant, Position<Frame>>{
      {lower_time, upper_time},
      {lower_degrees_of_freedom.position(),
       upper_degrees_of_freedom.position()},
      {lower_degrees_of_freedom.velocity(),
       upper_degrees_of_freedom.velocity()}};
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
    serialized_downsampling_parameters->set_max_dense_intervals(
        downsampling_parameters_->max_dense_intervals);
    downsampling_parameters_->tolerance.WriteToMessage(
        serialized_downsampling_parameters->mutable_tolerance());
  }
  message->set_number_of_dense_points(
      std::max<std::int64_t>(
          0, number_of_dense_points_ - number_of_points_to_skip_at_end));

  // Convert the |exact| vector into a set, and add the extremities.  This
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
  t.reserve(timeline_size);
  qx.reserve(timeline_size);
  qy.reserve(timeline_size);
  qz.reserve(timeline_size);
  px.reserve(timeline_size);
  py.reserve(timeline_size);
  pz.reserve(timeline_size);
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
    if (previous_instant.has_value()) {
      max_Δt = std::max(max_Δt, instant - *previous_instant);
    }
    previous_instant = instant;
  }

  // Times are exact.
  ZfpCompressor time_compressor(0);
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
}

}  // namespace internal_discrete_trajectory_segment
}  // namespace physics
}  // namespace principia
