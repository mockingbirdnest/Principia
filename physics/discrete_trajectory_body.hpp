
#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>
#include <vector>

#include "astronomy/epoch.hpp"
#include "base/array.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "numerics/fit_hermite_spline.hpp"
#include "quantities/si.hpp"
#include "zfp.h"

namespace principia {
namespace physics {
namespace internal_forkable {

using geometry::Instant;

template<typename Frame>
Instant const& ForkableTraits<DiscreteTrajectory<Frame>>::time(
    TimelineConstIterator const it) {
  return it->first;
}

template<typename Frame>
typename DiscreteTrajectoryIterator<Frame>::reference
DiscreteTrajectoryIterator<Frame>::operator*() const {
  auto const& it = this->current();
  return {it->first, it->second};
}

template<typename Frame>
std::optional<typename DiscreteTrajectoryIterator<Frame>::reference>
    DiscreteTrajectoryIterator<Frame>::operator->() const {
  auto const& it = this->current();
  return std::make_optional<reference>({it->first, it->second});
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame>*>
DiscreteTrajectoryIterator<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame> const*>
DiscreteTrajectoryIterator<Frame>::that() const {
  return this;
}

}  // namespace internal_forkable

namespace internal_discrete_trajectory {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using base::check_not_null;
using base::make_not_null_unique;
using base::not_null;
using base::UniqueArray;
using geometry::Displacement;
using numerics::FitHermiteSpline;
using quantities::si::Metre;
using quantities::si::Second;

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork()->time))
      << "NewForkWithCopy at nonexistent time " << time;

  auto const fork = this->NewFork(timeline_it);

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_.end()) {
    fork->timeline_.insert(++timeline_it, timeline_.end());
  }
  return fork;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithoutCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork()->time))
      << "NewForkWithoutCopy at nonexistent time " << time;

  return this->NewFork(timeline_it);
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkAtLast() {
  auto end = timeline_.end();
  if (timeline_.empty()) {
    return this->NewFork(end);
  } else {
    return this->NewFork(--end);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::AttachFork(
    not_null<std::unique_ptr<DiscreteTrajectory<Frame>>> fork) {
  CHECK(fork->is_root());
  CHECK(!this->Empty());

  auto& fork_timeline = fork->timeline_;
  auto const this_last = --this->end();

  // Determine if |fork| already has a point matching the end of this
  // trajectory.
  bool must_prepend;
  if (fork_timeline.empty()) {
    must_prepend = true;
  } else {
    CHECK_LE(this_last->time, fork_timeline.begin()->first);
    auto const it = fork_timeline.find(this_last->time);
    if (it == fork_timeline.end()) {
      must_prepend = true;
    } else {
      CHECK(it == fork_timeline.begin())
          << it->first << " " << this_last->time;
      must_prepend = false;
    }
  }

  // If needed, prepend to |fork| a copy of the last point of this trajectory.
  // This ensures that |fork| and this trajectory start and end, respectively,
  // with points at the same time (but possibly distinct degrees of freedom).
  if (must_prepend) {
    fork_timeline.emplace_hint(fork_timeline.begin(),
                               this_last->time,
                               this_last->degrees_of_freedom);
  }

  // Attach |fork| to this trajectory.
  this->AttachForkToCopiedBegin(std::move(fork));

  // Remove the first point of |fork| now that it is properly attached to its
  // parent: that point is either redundant (if it was prepended above) or wrong
  // (because we "trust" this trajectory more than |fork|).  The children that
  // might have been forked at the deleted point were relocated by
  // AttachForkToCopiedBegin.
  fork_timeline.erase(fork_timeline.begin());
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::DetachFork() {
  CHECK(!this->is_root());

  // Insert a new point in the timeline for the fork time.  It should go at the
  // beginning of the timeline.
  auto const fork_it = this->Fork();
  auto const begin_it = timeline_.emplace_hint(
      timeline_.begin(), fork_it->time, fork_it->degrees_of_freedom);
  CHECK(begin_it == timeline_.begin());

  // Detach this trajectory and tell the caller that it owns the pieces.
  return this->DetachForkWithCopiedBegin();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  CHECK(this->is_root() || time > this->Fork()->time)
       << "Append at " << time << " which is before fork time "
       << this->Fork()->time;

  if (!timeline_.empty() && timeline_.cbegin()->first == time) {
    LOG(WARNING) << "Append at existing time " << time
                 << ", time range = [" << this->front().time << ", "
                 << this->back().time << "]";
    return;
  }
  auto it = timeline_.emplace_hint(timeline_.end(),
                                   time,
                                   degrees_of_freedom);
  // Decrementing |end()| is much faster than incrementing |it|.  Don't ask.
  CHECK(--timeline_.end() == it)
      << "Append out of order at " << time << ", last time is "
      << (--timeline_.end())->first;
  if (downsampling_.has_value()) {
    if (timeline_.size() == 1) {
      downsampling_->SetStartOfDenseTimeline(timeline_.begin(), timeline_);
    } else {
      this->CheckNoForksBefore(this->back().time);
      downsampling_->increment_dense_intervals(timeline_);
      if (downsampling_->reached_max_dense_intervals()) {
        std::vector<TimelineConstIterator> dense_iterators;
        // This contains points, hence one more than intervals.
        dense_iterators.reserve(downsampling_->max_dense_intervals() + 1);
        for (TimelineConstIterator it =
                 downsampling_->start_of_dense_timeline();
             it != timeline_.end();
             ++it) {
          dense_iterators.push_back(it);
        }
        auto right_endpoints = FitHermiteSpline<Instant, Position<Frame>>(
            dense_iterators,
            [](auto&& it) -> auto&& { return it->first; },
            [](auto&& it) -> auto&& { return it->second.position(); },
            [](auto&& it) -> auto&& { return it->second.velocity(); },
            downsampling_->tolerance());
        if (right_endpoints.empty()) {
          right_endpoints.push_back(dense_iterators.end() - 1);
        }
        TimelineConstIterator left = downsampling_->start_of_dense_timeline();
        for (const auto& it_in_dense_iterators : right_endpoints) {
          TimelineConstIterator const right = *it_in_dense_iterators;
          timeline_.erase(++left, right);
          left = right;
        }
        downsampling_->SetStartOfDenseTimeline(left, timeline_);
      }
    }
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& time) {
  this->DeleteAllForksAfter(time);

  // Get an iterator denoting the first entry with time > |time|.  Remove that
  // entry and all the entries that follow it.  This preserves any entry with
  // time == |time|.
  auto const first_removed_in_timeline = timeline_.upper_bound(time);
  Instant const* const first_removed_time =
      first_removed_in_timeline == timeline_.end()
          ? nullptr
          : &first_removed_in_timeline->first;
  if (downsampling_.has_value()) {
    if (first_removed_time != nullptr &&
        *first_removed_time <= downsampling_->first_dense_time()) {
      // The start of the dense timeline will be invalidated.
      if (first_removed_in_timeline == timeline_.begin()) {
        // The timeline will be empty after erasing.
        downsampling_->SetStartOfDenseTimeline(timeline_.end(), timeline_);
      } else {
        // Further points will be appended to the last remaining point, so this
        // is where the dense timeline will begin.
        auto last_kept_in_timeline = first_removed_in_timeline;
        --last_kept_in_timeline;
        downsampling_->SetStartOfDenseTimeline(last_kept_in_timeline,
                                               timeline_);
      }
    }
  }
  timeline_.erase(first_removed_in_timeline, timeline_.end());
  if (downsampling_.has_value()) {
    downsampling_->RecountDenseIntervals(timeline_);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& time) {
  this->CheckNoForksBefore(time);

  // Get an iterator denoting the first entry with time >= |time|.  Remove all
  // the entries that precede it.  This preserves any entry with time == |time|.
  auto const first_kept_in_timeline = timeline_.lower_bound(time);
  if (downsampling_.has_value() &&
      (first_kept_in_timeline == timeline_.end() ||
       downsampling_->first_dense_time() < first_kept_in_timeline->first)) {
    // The start of the dense timeline will be invalidated.
    downsampling_->SetStartOfDenseTimeline(first_kept_in_timeline, timeline_);
  }
  timeline_.erase(timeline_.begin(), first_kept_in_timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::SetDownsampling(
    std::int64_t const max_dense_intervals,
    Length const& tolerance) {
  CHECK(this->is_root());
  CHECK(!downsampling_.has_value());
  downsampling_.emplace(
      max_dense_intervals, tolerance, timeline_.begin(), timeline_);
}
template<typename Frame>
void DiscreteTrajectory<Frame>::ClearDownsampling() {
  downsampling_.reset();
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_min() const {
  return this->Empty() ? InfiniteFuture : this->front().time;
}

template<typename Frame>
Instant DiscreteTrajectory<Frame>::t_max() const {
  return this->Empty() ? InfinitePast : this->back().time;
}

template<typename Frame>
Position<Frame> DiscreteTrajectory<Frame>::EvaluatePosition(
    Instant const& time) const {
  return GetInterpolation(time).Evaluate(time);
}

template<typename Frame>
Velocity<Frame> DiscreteTrajectory<Frame>::EvaluateVelocity(
    Instant const& time) const {;
  return GetInterpolation(time).EvaluateDerivative(time);
}

template<typename Frame>
DegreesOfFreedom<Frame> DiscreteTrajectory<Frame>::EvaluateDegreesOfFreedom(
    Instant const& time) const {
  auto const interpolation = GetInterpolation(time);
  return {interpolation.Evaluate(time), interpolation.EvaluateDerivative(time)};
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::DiscreteTrajectory*> const message,
    std::vector<DiscreteTrajectory<Frame>*> const& forks)
    const {
  CHECK(this->is_root());

  std::vector<DiscreteTrajectory<Frame>*> mutable_forks = forks;
  WriteSubTreeToMessage(message, mutable_forks);
  CHECK(std::all_of(mutable_forks.begin(),
                    mutable_forks.end(),
                    [](DiscreteTrajectory<Frame>* const fork) {
                      return fork == nullptr;
                    }));
}

template<typename Frame>
template<typename, typename>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& forks) {
  auto trajectory = make_not_null_unique<DiscreteTrajectory>();
  CHECK(std::all_of(forks.begin(),
                    forks.end(),
                    [](DiscreteTrajectory<Frame>** const fork) {
                      return fork != nullptr && *fork == nullptr;
                    }));
  trajectory->FillSubTreeFromMessage(message, forks);
  return trajectory;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*> DiscreteTrajectory<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame> const*>
DiscreteTrajectory<Frame>::that() const {
  return this;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_begin() const {
  return timeline_.begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_end() const {
  return timeline_.end();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_find(Instant const& time) const {
  return timeline_.find(time);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_lower_bound(Instant const& time) const {
  return timeline_.lower_bound(time);
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::timeline_size() const {
  return timeline_.size();
}

template<typename Frame>
DiscreteTrajectory<Frame>::Downsampling::Downsampling(
    std::int64_t const max_dense_intervals,
    Length const tolerance,
    TimelineConstIterator const start_of_dense_timeline,
    Timeline const& timeline)
    : max_dense_intervals_(max_dense_intervals),
      tolerance_(tolerance),
      start_of_dense_timeline_(start_of_dense_timeline) {
  RecountDenseIntervals(timeline);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::Downsampling::start_of_dense_timeline() const {
  return start_of_dense_timeline_;
}

template<typename Frame>
Instant const& DiscreteTrajectory<Frame>::Downsampling::first_dense_time()
    const {
  return start_of_dense_timeline_->first;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::SetStartOfDenseTimeline(
    TimelineConstIterator const value,
    Timeline const& timeline) {
  start_of_dense_timeline_ = value;
  RecountDenseIntervals(timeline);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::RecountDenseIntervals(
    Timeline const& timeline) {
  dense_intervals_ =
      std::distance(start_of_dense_timeline_, timeline.end()) - 1;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::increment_dense_intervals(
    Timeline const& timeline) {
  ++dense_intervals_;
  DCHECK_EQ(dense_intervals_,
            std::distance(start_of_dense_timeline_, timeline.end()) - 1);
}

template<typename Frame>
std::int64_t DiscreteTrajectory<Frame>::Downsampling::max_dense_intervals()
    const {
  return max_dense_intervals_;
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::Downsampling::reached_max_dense_intervals()
    const {
  return dense_intervals_ >= max_dense_intervals_;
}

template<typename Frame>
Length DiscreteTrajectory<Frame>::Downsampling::tolerance() const {
  return tolerance_;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Downsampling::WriteToMessage(
    not_null<serialization::DiscreteTrajectory::Downsampling*> message,
    Timeline const& timeline) const {
  if (start_of_dense_timeline_ == timeline.end()) {
    message->clear_start_of_dense_timeline();
  } else {
    first_dense_time().WriteToMessage(
        message->mutable_start_of_dense_timeline());
  }
  message->set_max_dense_intervals(max_dense_intervals_);
  tolerance_.WriteToMessage(message->mutable_tolerance());
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::Downsampling
DiscreteTrajectory<Frame>::Downsampling::ReadFromMessage(
    serialization::DiscreteTrajectory::Downsampling const& message,
    Timeline const& timeline) {
  TimelineConstIterator start_of_dense_timeline;
  if (message.has_start_of_dense_timeline()) {
    start_of_dense_timeline = timeline.find(
        Instant::ReadFromMessage(message.start_of_dense_timeline()));
    CHECK(start_of_dense_timeline != timeline.end());
  } else {
    start_of_dense_timeline = timeline.end();
  }
  return Downsampling(message.max_dense_intervals(),
                      Length::ReadFromMessage(message.tolerance()),
                      start_of_dense_timeline,
                      timeline);
}

inline void ZfpCompressDecompress(double const accuracy,
                                  int const input_size_in_bytes,
                                  zfp_field* const ifield,
                                  zfp_field* const ofield,
                                  int* const bytes,
                                  double* compression,
                                  int* const bpd,
                                  std::string* const out) {
  zfp_stream* zfp = zfp_stream_open(nullptr);
  if (accuracy == 0) {
    zfp_stream_set_reversible(zfp);
  } else {
    zfp_stream_set_accuracy(zfp, accuracy);
  }
  size_t bufsize = zfp_stream_maximum_size(zfp, ifield);
  void* buffer = malloc(bufsize);
  bitstream* stream = stream_open(buffer, bufsize);
  zfp_stream_set_bit_stream(zfp, stream);
  zfp_write_header(zfp, ifield, ZFP_HEADER_FULL);
  size_t izfpsize = zfp_compress(zfp, ifield);
  CHECK_LT(0, izfpsize);
  CHECK_EQ(izfpsize, zfp_stream_compressed_size(zfp));
  out->append(static_cast<char const*>(stream_data(stream)),
              stream_size(stream));

  zfp_stream_rewind(zfp);
  zfp_read_header(zfp, ofield, ZFP_HEADER_FULL);
  size_t ozfpsize = zfp_decompress(zfp, ofield);

  free(buffer);
  zfp_stream_close(zfp);

  *bytes = izfpsize;
  *compression = (double)input_size_in_bytes / izfpsize;
  *bpd = (8 * izfpsize) / (input_size_in_bytes / 8);
}

inline void ZfpWriteToMessage(double const accuracy,
                              const zfp_field* const field,
                              not_null<std::string*> const message) {
  std::unique_ptr<zfp_stream, std::function<void(zfp_stream*)>> const zfp(
      zfp_stream_open(/*stream=*/nullptr),
      [](zfp_stream* const zfp) { zfp_stream_close(zfp); });

  if (accuracy == 0) {
    zfp_stream_set_reversible(zfp.get());
  } else {
    zfp_stream_set_accuracy(zfp.get(), accuracy);
  }
  size_t const buffer_size = zfp_stream_maximum_size(zfp.get(), field);
  UniqueArray<std::uint8_t> const buffer(buffer_size);
  not_null<bitstream*> const stream =
      check_not_null(stream_open(buffer.data.get(), buffer_size));
  zfp_stream_set_bit_stream(zfp.get(), &*stream);

  zfp_write_header(zfp.get(), field, ZFP_HEADER_FULL);
  size_t const compressed_size = zfp_compress(zfp.get(), field);
  CHECK_LT(0, compressed_size);
  message->append(static_cast<char const*>(stream_data(stream)),
                  stream_size(stream));
}

inline void ZfpWriteToMessage2D(double const accuracy,
                                std::vector<double> const& v,
                                not_null<std::string*> const message) {
  constexpr int block = 4;
  auto const encoded = new double[(v.size() + block - 1) / block][block];
  for (int i = 0; i < (v.size() + block - 1) / block; ++i) {
    for (int j = 0; j < block; ++j) {
      if (block * i + j < v.size()) {
        encoded[i][j] = v[block * i + j];
      } else {
        // This will lead to poor compression at the end, but there is no
        // support for "ignored" data in zfp at this point.
        encoded[i][j] = 0;
      }
    }
  }

  // Beware nx and ny!  (And the Jabberwock, my son!)
  // See https://zfp.readthedocs.io/en/release0.5.5/tutorial.html#high-level-c-interface
  std::unique_ptr<zfp_field, std::function<void(zfp_field*)>> const field(
      zfp_field_2d(encoded,
                   /*type=*/zfp_type_double,
                   /*nx=*/block,
                   /*ny=*/(v.size() + block - 1) / block),
      [](zfp_field* const field) { zfp_field_free(field); });
  ZfpWriteToMessage(accuracy, field.get(), message);

  delete[] encoded;
}

inline void ZfpReadFromMessage(zfp_field* const field,
                               std::string_view& message) {
  std::unique_ptr<zfp_stream, std::function<void(zfp_stream*)>> const zfp(
      zfp_stream_open(/*stream=*/nullptr),
      [](zfp_stream* const zfp) { zfp_stream_close(zfp); });

  zfp_read_header(zfp.get(), field, ZFP_HEADER_FULL);
  size_t const field_size = zfp_field_size(field, /*size=*/nullptr);
  not_null<bitstream*> const stream = check_not_null(
      stream_open(const_cast<char*>(&message.front()), field_size));
  zfp_stream_set_bit_stream(zfp.get(), &*stream);

  size_t const uncompressed_size = zfp_decompress(zfp.get(), field);
  CHECK_LT(0, uncompressed_size);
  message.remove_prefix(field_size);
}

inline void ZfpReadFromMessage2D(int const size,
                                 std::vector<double>& v,
                                 std::string_view& message) {
  constexpr int block = 4;
  auto decoded = new double[(size + block - 1) / block][block];

  std::unique_ptr<zfp_field, std::function<void(zfp_field*)>> const field(
      zfp_field_2d(decoded,
                   /*type=*/zfp_type_double,
                   /*nx=*/block,
                   /*ny=*/(size + block - 1) / block),
      [](zfp_field* const field) { zfp_field_free(field); });

  ZfpReadFromMessage(field.get(), message);

  for (int i = 0; i < (size + block - 1) / block; ++i) {
    for (int j = 0; j < block; ++j) {
      if (block * i + j < size) {
        v[block * i + j] = decoded[i][j];
      } else {
        // Filler data.
      }
    }
  }
  delete[] decoded;
}

inline void ZpfExperiment2D(std::string_view const label,
                            double const accuracy,
                            std::vector<double> const& v,
                            std::string* const out) {
  auto input = new double[(v.size() + 3) / 4][4];
  auto output = new double[(v.size() + 3) / 4][4];
  for (int i = 0; i < (v.size() + 3) / 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (4 * i + j < v.size()) {
        input[i][j] = v[4 * i + j];
      } else {
        input[i][j] = 0;
      }
    }
  }

  zfp_type type = zfp_type_double;
  zfp_field* field = zfp_field_2d(input, type, 4, (v.size() + 3) / 4);
  zfp_field* ofield = zfp_field_2d(output, type, 4, (v.size() + 3) / 4);
  int bytes;
  double compression;
  int bpd;
  ZfpCompressDecompress(
      accuracy, 8 * v.size(), field, ofield, &bytes, &compression, &bpd, out);
  double max = 0;
  for (int i = 0; i < (v.size() + 3) / 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      max = std::max(max, std::abs(input[i][j] - output[i][j]));
    }
  }
  LOG(ERROR) << label << ": bytes: " << bytes
             << ": compression: " << compression << " bpd: " << bpd
             << " error: " << max;
  zfp_field_free(field);
  zfp_field_free(ofield);
  delete[] input;
  delete[] output;
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteSubTreeToMessage(
    not_null<serialization::DiscreteTrajectory*> const message,
    std::vector<DiscreteTrajectory<Frame>*>& forks) const {
  Forkable<DiscreteTrajectory, Iterator>::WriteSubTreeToMessage(message, forks);
  std::vector<double> t;
  std::vector<double> qx;
  std::vector<double> qy;
  std::vector<double> qz;
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  message->set_zfp_timeline_size(timeline_.size());
  std::string* const zfp_timeline = message->mutable_zfp_timeline();
  for (auto const& [instant, degrees_of_freedom] : timeline_) {
    auto const q = degrees_of_freedom.position() - Frame::origin;
    auto const p = degrees_of_freedom.velocity();
    t.push_back((instant - Instant{}) / Second);
    qx.push_back(q.coordinates().x / Metre);
    qy.push_back(q.coordinates().y / Metre);
    qz.push_back(q.coordinates().z / Metre);
    px.push_back(p.coordinates().x / (Metre / Second));
    py.push_back(p.coordinates().y / (Metre / Second));
    pz.push_back(p.coordinates().z / (Metre / Second));
  }
  //ZpfExperiment2D("t", 0, t, zfp_timeline);
  //ZpfExperiment2D("q.x", 10.0, qx, zfp_timeline);
  //ZpfExperiment2D("q.y", 10.0, qy, zfp_timeline);
  //ZpfExperiment2D("q.z", 10.0, qz, zfp_timeline);
  //ZpfExperiment2D("p.x", 0.1, px, zfp_timeline);
  //ZpfExperiment2D("p.y", 0.1, py, zfp_timeline);
  //ZpfExperiment2D("p.z", 0.1, pz, zfp_timeline);
  ZfpWriteToMessage2D(0, t, zfp_timeline);
  ZfpWriteToMessage2D(10.0, qx, zfp_timeline);
  ZfpWriteToMessage2D(10.0, qy, zfp_timeline);
  ZfpWriteToMessage2D(10.0, qz, zfp_timeline);
  ZfpWriteToMessage2D(0.1, px, zfp_timeline);
  ZfpWriteToMessage2D(0.1, py, zfp_timeline);
  ZfpWriteToMessage2D(0.1, pz, zfp_timeline);
  if (downsampling_.has_value()) {
    downsampling_->WriteToMessage(message->mutable_downsampling(), timeline_);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::FillSubTreeFromMessage(
    serialization::DiscreteTrajectory const& message,
    std::vector<DiscreteTrajectory<Frame>**> const& forks) {
  bool const is_pre_frobenius = message.timeline_size() > 0;
  if (is_pre_frobenius) {
    for (auto timeline_it = message.timeline().begin();
         timeline_it != message.timeline().end();
         ++timeline_it) {
      Append(Instant::ReadFromMessage(timeline_it->instant()),
             DegreesOfFreedom<Frame>::ReadFromMessage(
                 timeline_it->degrees_of_freedom()));
    }
  } else {
    int const size = message.zfp_timeline_size();
    std::vector<double> t(size);
    std::vector<double> qx(size);
    std::vector<double> qy(size);
    std::vector<double> qz(size);
    std::vector<double> px(size);
    std::vector<double> py(size);
    std::vector<double> pz(size);
    std::string_view zfp_timeline(message.zfp_timeline().data(),
                                  message.zfp_timeline().size());
    ZfpReadFromMessage2D(size, t, zfp_timeline);
    ZfpReadFromMessage2D(size, qx, zfp_timeline);
    ZfpReadFromMessage2D(size, qy, zfp_timeline);
    ZfpReadFromMessage2D(size, qz, zfp_timeline);
    ZfpReadFromMessage2D(size, px, zfp_timeline);
    ZfpReadFromMessage2D(size, py, zfp_timeline);
    ZfpReadFromMessage2D(size, pz, zfp_timeline);
    for (int i = 0; i < size; ++i) {
      Position<Frame> const q =
          Frame::origin +
          Displacement<Frame>({qx[i] * Metre, qy[i] * Metre, qz[i] * Metre});
      Velocity<Frame> const p({px[i] * (Metre / Second),
                               py[i] * (Metre / Second),
                               pz[i] * (Metre / Second)});
      Append(Instant() + t[i] * Second, DegreesOfFreedom<Frame>(q, p));
    }
  }
  if (message.has_downsampling()) {
    CHECK(this->is_root());
    downsampling_.emplace(
        Downsampling::ReadFromMessage(message.downsampling(), timeline_));
  }
  Forkable<DiscreteTrajectory, Iterator>::FillSubTreeFromMessage(message,
                                                                 forks);
}

template<typename Frame>
Hermite3<Instant, Position<Frame>> DiscreteTrajectory<Frame>::GetInterpolation(
    Instant const& time) const {
  CHECK_LE(t_min(), time);
  CHECK_GE(t_max(), time);
  // This is the upper bound of the interval upon which we will do the
  // interpolation.
  auto const upper = this->LowerBound(time);
  auto const lower = upper == this->begin() ? upper : --Iterator{upper};
  return Hermite3<Instant, Position<Frame>>{
      {lower->time, upper->time},
      {lower->degrees_of_freedom.position(),
       upper->degrees_of_freedom.position()},
      {lower->degrees_of_freedom.velocity(),
       upper->degrees_of_freedom.velocity()}};
}

}  // namespace internal_discrete_trajectory
}  // namespace physics
}  // namespace principia
