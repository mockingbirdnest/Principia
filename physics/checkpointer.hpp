
#pragma once

#include <functional>
#include <map>
#include <set>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/repeated_field.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

using base::not_null;
using geometry::Instant;
using quantities::Time;

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
// The |Message| must declare a nested message named |Checkpoint|, which must
// have a field named |time| of type |Point|.  There must be a repeated field of
// |Checkpoint|s in |Message|.
// This class is thread-safe.  The callbacks are not run under a lock.
template<typename Message>
class Checkpointer {
 public:
  // A function that writes an object to a checkpoint.  This function is
  // expected to capture the object being serialized.
  using Writer = std::function<void(not_null<typename Message::Checkpoint*>)>;

  // A function that reconstructs an object from a checkpoint as a side effect.
  // This function is expected to capture the object being deserialized.
  using Reader = std::function<Status(typename Message::Checkpoint const&)>;

  Checkpointer(Writer writer, Reader reader);

  // Returns the oldest checkpoint in this object, or +∞ if no checkpoint was
  // ever created.
  Instant oldest_checkpoint() const EXCLUDES(lock_);

  // Returns the newest checkpoint in this object, or -∞ if no checkpoint was
  // ever created.
  Instant newest_checkpoint() const EXCLUDES(lock_);

  // Returns the checkpoint at or immediately before |t|, or -∞ if no such
  // checkpoint exists.
  Instant checkpoint_at_or_before(Instant const& t) const EXCLUDES(lock_);

  // Returns all the checkpoints in this object.
  std::set<Instant> all_checkpoints() const EXCLUDES(lock_);

  // Returns all the checkpoints at or before |t|.
  std::set<Instant> all_checkpoints_at_or_before(Instant const& t) const
      EXCLUDES(lock_);

  // Creates a checkpoint at time |t|, which will be used to recreate the
  // timeline after |t|.  The checkpoint is constructed by calling the |Writer|
  // passed at construction.
  void WriteToCheckpoint(Instant const& t) EXCLUDES(lock_);

  // Same as above, but a checkpoint is only created if one was not created
  // recently, as specified by |max_time_between_checkpoints|.  Returns true iff
  // a new checkpoint was created.
  bool WriteToCheckpointIfNeeded(Instant const& t,
                                 Time const& max_time_between_checkpoints)
      EXCLUDES(lock_);

  // Calls the |Reader| passed at construction to reconstruct an object using
  // the oldest checkpoint.  Returns an error if this object contains no
  // checkpoint or if the |Reader| returns one.
  Status ReadFromOldestCheckpoint() const EXCLUDES(lock_);

  // Calls the |Reader| passed at construction to reconstruct an object using
  // the newest checkpoint.  Returns an error if this object contains no
  // checkpoint or if the |Reader| returns one.
  Status ReadFromNewestCheckpoint() const EXCLUDES(lock_);

  // Calls the |Reader| passed at construction to reconstruct an object using
  // the checkpoint at or immediately before |t|.  Returns an error if no such
  // checkpoint exists or if the |Reader| returns one.
  Status ReadFromCheckpointAtOrBefore(Instant const& t) const EXCLUDES(lock_);

  // Calls |reader| on the checkpoint at |t|.  Returns an error if there is no
  // such checkpoint or if |reader| returns one.
  Status ReadFromCheckpointAt(Instant const& t,
                              Reader const& reader) const EXCLUDES(lock_);

  // Calls |reader| on each of the checkpoints in this object, going backwards
  // from the most recent to the oldest.  Returns an error if |reader| returns
  // one.
  Status ReadFromAllCheckpointsBackwards(Reader const& reader) const
      EXCLUDES(lock_);

  void WriteToMessage(not_null<google::protobuf::RepeatedPtrField<
                          typename Message::Checkpoint>*> message) const
      EXCLUDES(lock_);
  static not_null<std::unique_ptr<Checkpointer>> ReadFromMessage(
      Writer writer,
      Reader reader,
      google::protobuf::RepeatedPtrField<typename Message::Checkpoint> const&
          message);

 private:
  void WriteToCheckpointLocked(Instant const& t)
      EXCLUSIVE_LOCKS_REQUIRED(lock_);

  mutable absl::Mutex lock_;
  Writer const writer_;
  Reader const reader_;

  // The time field of the Checkpoint message may or may not be set.  The map
  // key is the source of truth.
  std::map<Instant, typename Message::Checkpoint> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
