#pragma once

#include "physics/checkpointer.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

template<typename Message, typenameTypes...>
Checkpointer<Message, Types...>::Checkpointer(Filler filler,
                                        Reader reader,
                                        Writer writer)
    : filler_(std::move(filler)),
      reader_(std::move(reader)),
      writer_(std::move(writer)) {}

template<typename Message, typename... Types>
void Checkpointer<Message, ... Types>::CreateIfNeeded(
    Instant const& t,
    Time const& max_time_between_checkpoints) {}

template<typename Message, typename... Types>
void Checkpointer<Message, ... Types>::CreateUnconditionally(Instant const& t) {
}

template<typename Message, typenameTypes...>
bool Checkpointer<Message, Types...>::HasBefore(Instant const& t) {
  return false;
}

template<typename Message, typenameTypes...>
void Checkpointer<Message, Types...>::ForgetBefore(Instant const& t) {}

template<typename Message, typename... Types>
void Checkpointer<Message, ... Types>::WriteToMessage(
    not_null<Message*> message) {}

template<typename Message, typename... Types>
void Checkpointer<Message, ... Types>::ReadFromMessage(Message const& message) {
}

}  // namespace internal_checkpointer
}  // namespace physics
}  // namespace principia
