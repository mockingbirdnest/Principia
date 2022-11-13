#pragma once

#include "physics/checkpointer.hpp"

#include <map>
#include <vector>
#include <utility>

#include "absl/container/btree_set.h"
#include "base/status_utilities.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

using geometry::InfiniteFuture;
using geometry::InfinitePast;

template<typename Message>
Checkpointer<Message>::Checkpointer(Writer writer, Reader reader)
    : writer_(std::move(writer)),
      reader_(std::move(reader)) {}

template<typename Message>
Instant Checkpointer<Message>::oldest_checkpoint() const {
  absl::ReaderMutexLock l(&lock_);
  if (checkpoints_.empty()) {
    return InfiniteFuture;
  }
  return checkpoints_.cbegin()->first;
}

template<typename Message>
Instant Checkpointer<Message>::newest_checkpoint() const {
  absl::ReaderMutexLock l(&lock_);
  if (checkpoints_.empty()) {
    return InfinitePast;
  }
  return checkpoints_.crbegin()->first;
}

template<typename Message>
Instant Checkpointer<Message>::checkpoint_at_or_after(Instant const& t) const {
  absl::ReaderMutexLock l(&lock_);
  // |it| denotes an entry equal to or greater than |t| (or end).
  auto const it = checkpoints_.lower_bound(t);
  if (it == checkpoints_.cend()) {
    return InfiniteFuture;
  }
  return it->first;
}

template<typename Message>
Instant Checkpointer<Message>::checkpoint_at_or_before(Instant const& t) const {
  absl::ReaderMutexLock l(&lock_);
  // |it| denotes an entry strictly greater than |t| (or end).
  auto const it = checkpoints_.upper_bound(t);
  if (it == checkpoints_.cbegin()) {
    return InfinitePast;
  }
  return std::prev(it)->first;
}

template<typename Message>
absl::btree_set<Instant> Checkpointer<Message>::all_checkpoints() const {
  absl::ReaderMutexLock l(&lock_);
  absl::btree_set<Instant> result;
  std::transform(
      checkpoints_.cbegin(),
      checkpoints_.cend(),
      std::inserter(result, result.end()),
      [](typename CheckpointsByTime::value_type const& pair) {
        return pair.first;
      });
  return result;
}

template<typename Message>
absl::btree_set<Instant> Checkpointer<Message>::all_checkpoints_at_or_before(
    Instant const& t) const {
  absl::ReaderMutexLock l(&lock_);
  // |it| denotes an entry strictly greater than |t| (or end).
  auto const it = checkpoints_.upper_bound(t);
  absl::btree_set<Instant> result;
  std::transform(
      checkpoints_.cbegin(),
      it,
      std::inserter(result, result.end()),
      [](typename CheckpointsByTime::value_type const& pair) {
        return pair.first;
      });
  return result;
}

template<typename Message>
absl::btree_set<Instant> Checkpointer<Message>::all_checkpoints_between(
    Instant const& t1,
    Instant const& t2) const {
  if (t2 < t1) {
    return absl::btree_set<Instant>();
  }

  absl::ReaderMutexLock l(&lock_);
  // |it1| denotes an entry greater or equal to |t1| (or end).
  auto const it1 = checkpoints_.lower_bound(t1);
  // |it2| denotes an entry strictly greater than |t2| (or end).
  auto const it2 = checkpoints_.upper_bound(t2);
  absl::btree_set<Instant> result;
  std::transform(
      it1,
      it2,
      std::inserter(result, result.end()),
      [](typename CheckpointsByTime::value_type const& pair) {
        return pair.first;
      });
  return result;
}

template<typename Message>
void Checkpointer<Message>::WriteToCheckpoint(Instant const& t) {
  absl::MutexLock l(&lock_);
  WriteToCheckpointLocked(t);
}

template<typename Message>
bool Checkpointer<Message>::WriteToCheckpointIfNeeded(
    Instant const& t,
    Time const& max_time_between_checkpoints) {
  absl::MutexLock l(&lock_);
  if (checkpoints_.empty() ||
      max_time_between_checkpoints < t - checkpoints_.crbegin()->first) {
    WriteToCheckpointLocked(t);
    return true;
  }
  return false;
}

template<typename Message>
absl::Status Checkpointer<Message>::ReadFromOldestCheckpoint() const {
  typename Message::Checkpoint const* checkpoint = nullptr;
  {
    absl::ReaderMutexLock l(&lock_);
    if (checkpoints_.empty()) {
      return absl::NotFoundError("No checkpoint");
    }
    checkpoint = &checkpoints_.cbegin()->second;
  }
  return reader_(*checkpoint);
}

template<typename Message>
absl::Status Checkpointer<Message>::ReadFromNewestCheckpoint() const {
  typename Message::Checkpoint const* checkpoint = nullptr;
  {
    absl::ReaderMutexLock l(&lock_);
    if (checkpoints_.empty()) {
      return absl::NotFoundError("No checkpoint");
    }
    checkpoint = &checkpoints_.crbegin()->second;
  }
  return reader_(*checkpoint);
}

template<typename Message>
absl::Status Checkpointer<Message>::ReadFromCheckpointAtOrBefore(
    Instant const& t) const {
  typename Message::Checkpoint const* checkpoint = nullptr;
  {
    absl::ReaderMutexLock l(&lock_);
    // |it| denotes an entry strictly greater than |t| (or end).
    auto const it = checkpoints_.upper_bound(t);
    if (it == checkpoints_.cbegin()) {
      return absl::NotFoundError("No checkpoint");
    }
    checkpoint = &std::prev(it)->second;
  }
  return reader_(*checkpoint);
}

template<typename Message>
absl::Status Checkpointer<Message>::ReadFromCheckpointAt(
    Instant const& t,
    Reader const& reader) const {
  typename CheckpointsByTime::const_iterator it;
  {
    absl::ReaderMutexLock l(&lock_);
    it = checkpoints_.find(t);
    if (it == checkpoints_.end()) {
      return absl::NotFoundError("No checkpoint found");
    }
  }
  return reader(it->second);
}

template<typename Message>
absl::Status Checkpointer<Message>::ReadFromCheckpointAt(
    Instant const& t) const {
  return ReadFromCheckpointAt(t, reader_);
}

template<typename Message>
void Checkpointer<Message>::WriteToMessage(
    not_null<google::protobuf::RepeatedPtrField<typename Message::Checkpoint>*>
        message) const {
  absl::ReaderMutexLock l(&lock_);
  for (const auto [time, checkpoint] : checkpoints_) {
    typename Message::Checkpoint* const message_checkpoint = message->Add();
    *message_checkpoint = checkpoint;
    time.WriteToMessage(message_checkpoint->mutable_time());
  }
}

template<typename Message>
not_null<std::unique_ptr<Checkpointer<Message>>>
Checkpointer<Message>::ReadFromMessage(
    Writer writer,
    Reader reader,
    google::protobuf::RepeatedPtrField<typename Message::Checkpoint> const&
        message) {
  auto checkpointer =
      std::make_unique<Checkpointer>(std::move(writer), std::move(reader));
  for (const auto& checkpoint : message) {
    Instant const time = Instant::ReadFromMessage(checkpoint.time());
    checkpointer->checkpoints_.emplace(time, checkpoint);
  }
  return std::move(checkpointer);
}

template<typename Message>
void Checkpointer<Message>::WriteToCheckpointLocked(Instant const& t) {
  lock_.AssertHeld();
  CHECK(!checkpoints_.contains(t)) << t;
  auto const it = checkpoints_.emplace_hint(
      checkpoints_.end(), t, typename Message::Checkpoint());
  lock_.Unlock();
  writer_(&it->second);
  lock_.Lock();
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
