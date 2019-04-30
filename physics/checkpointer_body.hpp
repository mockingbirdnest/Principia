
#pragma once

#include "physics/checkpointer.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

template<typename Object, typename Message>
Checkpointer<Object, Message>::Checkpointer(Reader reader, Writer writer)
    : reader_(std::move(reader)),
      writer_(std::move(writer)) {}

template<typename Object, typename Message>
void Checkpointer<Object, Message>::CreateIfNeeded(
    Instant const& t,
    Time const& max_time_between_checkpoints) {
  if (checkpoints_.empty() ||
      max_time_between_checkpoints < t - checkpoints_.crbegin()->first) {
    CreateUnconditionally(t);
  }
}

template<typename Object, typename Message>
void Checkpointer<Object, Message>::CreateUnconditionally(Instant const& t) {
  auto const it = checkpoints_.emplace_hint(checkpoints_.end(), t, Message);
  writer_(&it->second);
}

template<typename Object, typename Message>
Instant const& Checkpointer<Object, Message>::OldestCheckpointTime() const {
  if (checkpoints_.empty()) {
    return InfinitePast;
  } else {
    return checkpoints_.begin()->first;
  }
}

template<typename Object, typename Message>
void Checkpointer<Object, Message>::ForgetBefore(Instant const& t) {
  auto const it = checkpoints_.upper_bound(t);
  CHECK(it == checkpoints_.end() || t < it->first);
  checkpoints_.erase(checkpoints_.begin(), it);
}

template<typename Object, typename Message>
void Checkpointer<Object, Message>::WriteToMessage(
    not_null<Message*> const message) {
  writer_(message);
}

template<typename Object, typename Message>
void Checkpointer<Object, Message>::ReadFromMessage(
    Message const& message,
    not_null<Object*> const object) {
  reader_(message, *object);
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
