
#pragma once

#include "physics/checkpointer.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

template<typename Message>
Checkpointer<Message>::Checkpointer(Reader reader, Writer writer)
    : reader_(std::move(reader)),
      writer_(std::move(writer)) {}

template<typename Message>
bool Checkpointer<Message>::CreateIfNeeded(
    Instant const& t,
    Time const& max_time_between_checkpoints) {
  if (checkpoints_.empty() ||
      max_time_between_checkpoints < t - checkpoints_.crbegin()->first) {
    CreateUnconditionally(t);
    return true;
  }
  return false;
}

template<typename Message>
void Checkpointer<Message>::CreateUnconditionally(Instant const& t) {
  auto const it = checkpoints_.emplace_hint(checkpoints_.end(), t, Message());
  writer_(&it->second);
}

template<typename Message>
void Checkpointer<Message>::ForgetBefore(Instant const& t) {
  auto const it = checkpoints_.upper_bound(t);
  CHECK(it == checkpoints_.end() || t < it->first);
  checkpoints_.erase(checkpoints_.begin(), it);
}

template<typename Message>
Instant Checkpointer<Message>::WriteToMessage(
    not_null<Message*> const message) const {
  if (checkpoints_.empty()) {
    // TODO(phl): declare this next to Instant.
    static Instant infinite_future = Instant() + quantities::Infinity<Time>();
    return infinite_future;
  } else {
    message->MergeFrom(checkpoints_.cbegin()->second);
    return checkpoints_.cbegin()->first;
  }
}

template<typename Message>
void Checkpointer<Message>::ReadFromMessage(Instant const& t,
                                            Message const& message) {
  checkpoints_.clear();
  if (reader_(message)) {
    CreateUnconditionally(t);
  }
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
