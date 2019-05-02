
#pragma once

#include <functional>
#include <tuple>

#include "absl/synchronization/mutex.h"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
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
// the timeline may be reconstructed as needed based on the checkpoint.
// Checkpoints must be created at regular intervals because they are dropped by
// ForgetBefore: this ensures that there is always a sufficient old checkpoint
// available the next time serialization is performed.
// Logically checkpoints would serialize to/deserialize from a specific message,
// but for historical reasons they just fill some fields of a message that may
// contain other information.
// This class is thread-safe.
template<typename Message>
class Checkpointer {
 public:
  // A function that reconstructs an object from a checkpoint as a side effect.
  // It must return true iff the given message actually contained a checkpoint.
  // If it returns false, it must leave the object in a state corresponding to
  // its default initialization.  This function is expected to capture the
  // object being deserialized.
  using Reader = std::function<bool(Message const&)>;

  // A function that writes an object to a checkpoint.  This function is
  // expected to capture the object being serialized.
  using Writer = std::function<void(not_null<Message*>)>;

  Checkpointer(Reader reader, Writer writer);

  // Creates a checkpoint at time |t|, which will be used to recreate the
  // timeline after |t|.  The checkpoint is constructed by calling the |Writer|
  // passed at construction.
  void CreateUnconditionally(Instant const& t) EXCLUDES(lock_);

  // Same as above, but a checkpoint is only created if one was not created
  // recently, as specified by |max_time_between_checkpoints|.
  bool CreateIfNeeded(Instant const& t,
                      Time const& max_time_between_checkpoints) EXCLUDES(lock_);

  // Removes all checkpoints for times strictly less than |t|.
  void ForgetBefore(Instant const& t) EXCLUDES(lock_);

  // If there exist a checkpoint, writes the oldest checkpoint to the |message|
  // using protocol buffer merging and returns its time.  Otherwise returns +∞.
  // The time returned by this function should be serialized and passed to
  // |ReadFromMessage| when deserializing to ensure that checkpoints are
  // preserved across serialization/deserialization cycles.
  Instant WriteToMessage(not_null<Message*> message) const EXCLUDES(lock_);

  // Clears all the checkpoints in this checkpointer, and calls the |Reader|
  // passed at construction to reconstruct the object from |message|.  If the
  // |Reader| returns true (i.e., there was a checkpoint in the |message|),
  // creates a new checkpoint in this checkpointer at the given time.
  void ReadFromMessage(Instant const& t,
                       Message const& message) EXCLUDES(lock_);

 private:
  void CreateUnconditionallyLocked(Instant const& t)
      EXCLUSIVE_LOCKS_REQUIRED(lock_);

  mutable absl::Mutex lock_;
  Reader const reader_;
  Writer const writer_;
  std::map<Instant, Message> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
