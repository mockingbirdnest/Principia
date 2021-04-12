
#pragma once

#include "physics/checkpointer.hpp"

#include <vector>

#include "astronomy/epoch.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

using astronomy::InfiniteFuture;
using astronomy::InfinitePast;
using base::Error;

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
std::set<Instant> Checkpointer<Message>::all_checkpoints() const {
  absl::ReaderMutexLock l(&lock_);
  std::set<Instant> result;
  std::transform(
      checkpoints_.cbegin(),
      checkpoints_.cend(),
      std::inserter(result, result.end()),
      [](std::pair<Instant, typename Message::Checkpoint> const& pair) {
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
Status Checkpointer<Message>::ReadFromOldestCheckpoint() const {
  typename Message::Checkpoint const* checkpoint = nullptr;
  {
    absl::ReaderMutexLock l(&lock_);
    if (checkpoints_.empty()) {
      return Status(Error::NOT_FOUND, "No checkpoint");
    }
    checkpoint = &checkpoints_.cbegin()->second;
  }
  return reader_(*checkpoint);
}

template<typename Message>
Status Checkpointer<Message>::ReadFromNewestCheckpoint() const {
  typename Message::Checkpoint const* checkpoint = nullptr;
  {
    absl::ReaderMutexLock l(&lock_);
    if (checkpoints_.empty()) {
      return Status(Error::NOT_FOUND, "No checkpoint");
    }
    checkpoint = &checkpoints_.crbegin()->second;
  }
  return reader_(*checkpoint);
}

template<typename Message>
Status Checkpointer<Message>::ReadFromCheckpointAt(Instant const& checkpoint,
                                                   Reader const& reader) const {
  typename std::map<Instant, Message::Checkpoint>::const_iterator it;
  {
    absl::ReaderMutexLock l(&lock_);
    it = checkpoints_.find(checkpoint);
    if (it == checkpoints_.end()) {
      return Status(Error::NOT_FOUND, "No checkpoint found");
    }
  }
  return reader(it->second);
}

template<typename Message>
Status Checkpointer<Message>::ReadFromAllCheckpointsBackwards(
    Reader const& reader) const {
  // We'll be running the callback without the lock, so we take a snapshot of
  // the checkpoints.
  std::vector<not_null<typename Message::Checkpoint const*>> checkpoints;
  {
    absl::ReaderMutexLock l(&lock_);
    for (auto it = checkpoints_.crbegin(); it != checkpoints_.crend(); ++it) {
      checkpoints.push_back(&it->second);
    }
  }
  for (auto const checkpoint : checkpoints) {
    RETURN_IF_ERROR(reader(*checkpoint));
  }
  return Status::OK;
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
  auto const it = checkpoints_.emplace_hint(
      checkpoints_.end(), t, typename Message::Checkpoint());
  lock_.Unlock();
  writer_(&it->second);
  lock_.Lock();
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
