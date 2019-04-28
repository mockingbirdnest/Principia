
#pragma once

#include <functional>
#include <tuple>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp

namespace principia {
namespace physics {
namespace internal_checkpointer {

using base::not_null;
using geometry::Instant;
using quantities::Time;

template<typename Message, typename... Types>
class Checkpointer {
 public:
  using Checkpoint = std::tuple<Types...>;
  using Filler = std::function<Checkpoint(T const&)>;
  using Reader = std::function<Checkpoint(Message const&)>;
  using Writer = std::function<void(Checkpoint const&, Message*)>;

  Checkpointer(Filler filler,
               Reader reader,
               Writer writer);

  void CreateIfNeeded(Instant const& t,
                      Time const& max_time_between_checkpoints);
  void CreateUnconditionally(Instant const& t);

  Instant const& OldestCheckpoint();
  bool HasBefore(Instant const& t);
  void ForgetBefore(Instant const& t);

  void WriteToMessage(not_null<Message*> message);
  static void ReadFromMessage(Message const& message);

 private:
  Filler const filler_;
  Reader const reader_;
  Writer const writer_;
  std::map<Instant, Checkpoint> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
