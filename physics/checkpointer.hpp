
#pragma once

#include <functional>
#include <tuple>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

using base::not_null;
using geometry::Instant;
using quantities::Time;

//TODO(phl):comment.  Thread safety?
template<typename Message>
class Checkpointer {
 public:
  using Reader = std::function<void(Message const&)>;
  using Writer = std::function<void(not_null<Message*>)>;

  Checkpointer(Reader reader, Writer writer);

  void CreateIfNeeded(Instant const& t,
                      Time const& max_time_between_checkpoints);
  void CreateUnconditionally(Instant const& t);

  void ForgetBefore(Instant const& t);

  Instant WriteToMessage(not_null<Message*> message) const;
  void ReadFromMessage(Instant const& t,
                       Message const& message);

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
