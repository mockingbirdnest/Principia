
#pragma once

#include "physics/checkpointer.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

template<typename T, typename Message>
Checkpointer<T, Message>::Checkpointer(Reader reader, Writer writer)
    : reader_(std::move(reader)),
      writer_(std::move(writer)) {}

template<typename T, typename Message>
void Checkpointer<T, Message>::CreateIfNeeded(
    Instant const& t,
    Time const& max_time_between_checkpoints) {
  if (checkpoints_.empty() ||
      max_time_between_checkpoints < t - checkpoints_.crbegin()->first) {
    CreateUnconditionally(t);
  }
}

template<typename T, typename Message>
void Checkpointer<T, Message>::CreateUnconditionally(Instant const& t) {
  auto const it = checkpoints_.emplace_hint(checkpoints_.end(), t, Message);
  writer_(&it->second);
}

template<typename T, typename Message>
Instant const& Checkpointer<T, Message>::OldestCheckpointTime() const {
  if (checkpoints_.empty()) {
    return InfinitePast;
  } else {
    return checkpoints_.begin()->first;
  }
}

template<typename T, typename Message>
void Checkpointer<T, Message>::ForgetBefore(Instant const& t) {
  auto const it = checkpoints_.upper_bound(t);
  CHECK(it == checkpoints_.end() || t < it->first);
  checkpoints_.erase(checkpoints_.begin(), it);
}

template<typename T, typename Message>
void Checkpointer<T, Message>::WriteToMessage(
    not_null<Message*> message) {}

template<typename T, typename Message>
void Checkpointer<T, Message>::ReadFromMessage(Message const& message) {
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
