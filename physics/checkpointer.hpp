
#pragma once

#include <functional>
#include <map>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
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
// checkpoint, followed by the checkpoint itself.  When deserializing, the
// timeline may be reconstructed as needed based on the checkpoint.
// Checkpoints must be created at regular intervals to ensure that the timeline
// may be reconstructed fast enough.
// Logically checkpoints would serialize to/deserialize from a specific message,
// but for historical reasons they just fill some fields of a message that may
// contain other information.
//TODO(phl):Fix comments
// This class is thread-safe.
template<typename Message>
class Checkpointer {
 public:
  // A function that reconstructs an object from a checkpoint as a side effect.
  // It must return true iff the given message actually contained a checkpoint.
  // If it returns false, it must leave the object in a state corresponding to
  // its default initialization.  This function is expected to capture the
  // object being deserialized.
  using Reader = std::function<void(typename Message::Checkpoint const&)>;

  // A function that writes an object to a checkpoint.  This function is
  // expected to capture the object being serialized.
  using Writer = std::function<void(not_null<typename Message::Checkpoint*>)>;

  Checkpointer(Reader reader,
               Writer writer);

  Instant oldest_checkpoint() const;

  // Creates a checkpoint at time |t|, which will be used to recreate the
  // timeline after |t|.  The checkpoint is constructed by calling the |Writer|
  // passed at construction.
  void CreateUnconditionally(Instant const& t) EXCLUDES(lock_);

  // Same as above, but a checkpoint is only created if one was not created
  // recently, as specified by |max_time_between_checkpoints|.
  bool CreateIfNeeded(Instant const& t,
                      Time const& max_time_between_checkpoints) EXCLUDES(lock_);

  // If there exists a checkpoint, writes the oldest checkpoint to the |message|
  // using protocol buffer merging and returns its time.  Otherwise returns +∞.
  // The time returned by this function should be serialized and passed to
  // |ReadFromMessage| when deserializing to ensure that checkpoints are
  // preserved across serialization/deserialization cycles.
  void WriteToMessage(not_null<google::protobuf::RepeatedPtrField<
                          typename Message::Checkpoint>*> message) const
      EXCLUDES(lock_);

  // Clears all the checkpoints in this checkpointer, and calls the |Reader|
  // passed at construction to reconstruct the object from |message|.  If the
  // |Reader| returns true (i.e., there was a checkpoint in the |message|),
  // creates a new checkpoint in this checkpointer at the given time.
  static not_null<std::unique_ptr<Checkpointer>> ReadFromMessage(
      Reader reader,
      Writer writer,
      google::protobuf::RepeatedPtrField<typename Message::Checkpoint> const&
          message);

 private:
  void CreateUnconditionallyLocked(Instant const& t)
      EXCLUSIVE_LOCKS_REQUIRED(lock_);

  mutable absl::Mutex lock_;
  Reader const reader_;
  Writer const writer_;

  // The time field of the Checkpoint message may or may not be set.  The map
  // key is the source of truth.
  std::map<Instant, typename Message::Checkpoint> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
