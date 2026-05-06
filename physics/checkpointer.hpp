#pragma once

#include <functional>
#include <memory>

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/status/status.h"
#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "google/protobuf/arena.h"
#include "google/protobuf/repeated_field.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace _checkpointer {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::quantities::_quantities;

// The checkpointer helps with compact serialization of timelines, i.e., classes
// that associate some data with distinct instants.  The naïve implementation of
// serialization for timelines would write all the pairs (time, data) but that
// would potentially result in large saves that would be expensive to read and
// write.
// Instead, this class creates checkpoints that encapsulate all information
// needed to reconstruct the timeline after a give point in time.  When
// serializing a timeline, the pairs (time, data) are written up to the oldest
// checkpoint, followed by the checkpoints themselves.  When deserializing, the
// timeline may be reconstructed as needed based on the checkpoints.
// Checkpoints must be created at regular intervals to ensure that the timeline
// may be reconstructed fast enough.
// The `Message` must declare a nested message named `Checkpoint`, which must
// have a field named `time` of type `Point`.  There must be a repeated field of
// `Checkpoint`s in `Message`.
// This class is thread-safe.  The callbacks are not run under a lock.
template<typename Message>
class Checkpointer {
 public:
  // A function that writes an object to a checkpoint.  This function is
  // expected to capture the object being serialized.
  using Writer = std::function<void(not_null<typename Message::Checkpoint*>)>;

  // A function that reconstructs an object from a checkpoint as a side effect.
  // This function is expected to capture the object being deserialized.
  using Reader =
      std::function<absl::Status(typename Message::Checkpoint const&)>;

// Saves produced by [Лефшец, Leibniz[ had a ton of empty segments at the
// beginning of the non-collapsible part of their checkpoints, causing storage
// to explode: we saw a vessel with 1.8 million checkpoints, about 10'000 of
// which were non-empty. To fix these saves, we rewrite the checkpoints.  That's
// what we get for writing quadratic code.
  using Rewriter = std::function<typename Message::Checkpoint(
      typename Message::Checkpoint const&)>;

  Checkpointer(Writer writer, Reader reader);

  // Updates the writer of this object.
  void set_writer(Writer writer) EXCLUDES(lock_);

  // The number of checkpoints held by this object.  Only use for debugging and
  // logging.
  std::int64_t size() const EXCLUDES(lock_);

  // Returns the oldest checkpoint in this object, or +∞ if no checkpoint was
  // ever created.
  Instant oldest_checkpoint() const EXCLUDES(lock_);

  // Returns the newest checkpoint in this object, or -∞ if no checkpoint was
  // ever created.
  Instant newest_checkpoint() const EXCLUDES(lock_);

  // Returns the checkpoint at or immediately after `t`, or +∞ if no such
  // checkpoint exists.
  Instant checkpoint_at_or_after(Instant const& t) const EXCLUDES(lock_);

  // Returns the checkpoint at or immediately before `t`, or -∞ if no such
  // checkpoint exists.
  Instant checkpoint_at_or_before(Instant const& t) const EXCLUDES(lock_);

  // Returns all the checkpoints in this object.
  absl::btree_set<Instant> all_checkpoints() const EXCLUDES(lock_);

  // Returns all the checkpoints at or before `t`.
  absl::btree_set<Instant> all_checkpoints_at_or_before(Instant const& t) const
      EXCLUDES(lock_);

  // Returns all the checkpoints in interval [t1, t2].
  absl::btree_set<Instant> all_checkpoints_between(Instant const& t1,
                                                   Instant const& t2) const
      EXCLUDES(lock_);

  // Creates a checkpoint at time `t`, which will be used to recreate the
  // timeline after `t`.  The checkpoint is constructed by calling the `Writer`
  // passed at construction.
  void WriteToCheckpoint(Instant const& t) EXCLUDES(lock_);

  // Same as above, but a checkpoint is only created if one was not created
  // recently, as specified by `max_time_between_checkpoints`.  Returns true iff
  // a new checkpoint was created.
  bool WriteToCheckpointIfNeeded(Instant const& t,
                                 Time const& max_time_between_checkpoints)
      EXCLUDES(lock_);

  // Calls the `Reader` passed at construction to reconstruct an object using
  // the oldest checkpoint.  Returns an error if this object contains no
  // checkpoint or if the `Reader` returns one.
  absl::Status ReadFromOldestCheckpoint() const EXCLUDES(lock_);

  // Calls the `Reader` passed at construction to reconstruct an object using
  // the newest checkpoint.  Returns an error if this object contains no
  // checkpoint or if the `Reader` returns one.
  absl::Status ReadFromNewestCheckpoint() const EXCLUDES(lock_);

  // Calls the `Reader` passed at construction to reconstruct an object using
  // the checkpoint at or immediately before `t`.  Returns an error if no such
  // checkpoint exists or if the `Reader` returns one.
  absl::Status ReadFromCheckpointAtOrBefore(Instant const& t) const
      EXCLUDES(lock_);

  // Calls `reader` on the checkpoint at `t`.  Returns an error if there is no
  // such checkpoint or if `reader` returns one.
  absl::Status ReadFromCheckpointAt(Instant const& t,
                                    Reader const& reader) const EXCLUDES(lock_);

  // Same as above, but uses the reader passed at construction.
  absl::Status ReadFromCheckpointAt(Instant const& t) const EXCLUDES(lock_);

  void WriteToMessage(not_null<google::protobuf::RepeatedPtrField<
                          typename Message::Checkpoint>*> message) const
      EXCLUDES(lock_);
  // `rewriter` may be null, in which case checkpoints are used as-is.
  static not_null<std::unique_ptr<Checkpointer>> ReadFromMessage(
      Writer writer,
      Reader reader,
      Rewriter rewriter,
      google::protobuf::RepeatedPtrField<typename Message::Checkpoint> const&
          message);

 private:
  using CheckpointsByTime = absl::btree_map<
      Instant,
      google::protobuf::Arena::UniquePtr<typename Message::Checkpoint>>;

  void WriteToCheckpointLocked(Instant const& t)
      ABSL_EXCLUSIVE_LOCKS_REQUIRED(lock_);

  mutable absl::Mutex lock_;
  Writer writer_;
  Reader const reader_;

  // The time field of the Checkpoint message may or may not be set.  The map
  // key is the source of truth.
  CheckpointsByTime checkpoints_;

  // The checkpoints held in the `checkpoints_` map are all allocated in this
  // arena.
  google::protobuf::Arena arena_;
};

}  // namespace internal

using internal::Checkpointer;

}  // namespace _checkpointer
}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
