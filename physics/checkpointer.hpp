
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

template<typename Object, typename Message>
class Checkpointer {
 public:
  using Reader = std::function<void(Message const&, Object&)>;
  using Writer = std::function<void(not_null<Message>*)>;

  Checkpointer(Reader reader, Writer writer);

  void CreateIfNeeded(Instant const& t,
                      Time const& max_time_between_checkpoints);
  void CreateUnconditionally(Instant const& t);

  Instant const& OldestCheckpointTime() const;
  void ForgetBefore(Instant const& t);

  void WriteToMessage(not_null<Message*> message);
  static void ReadFromMessage(Message const& message,
                              not_null<Object*> object);

 private:
  Reader const reader_;
  Writer const writer_;
  std::map<Instant, Message> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia

#include "physics/checkpointer_body.hpp"
