#pragma once

#include <functional>
#include <tuple>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_checkpointer {

using base::not_null;
using geometry::Instant;

template<typename T, typename... Types>
class Checkpointer {
 public:
  using Checkpoint = std::tuple<Types...>;
  using Filler = std::function<Checkpoint(T const&)>;

  Checkpointer(not_null<T const*> t, Filler filler);

  Checkpoint const& Get();
  void GetIfNeeded();

  bool HasBefore(Instant const& t);

  void ForgetBefore(Instant const& t);

  void WriteToMessage();
  static void ReadFromMessage();

private:
  std::vector<Checkpoint> checkpoints_;
};

}  // namespace internal_checkpointer

using internal_checkpointer::Checkpointer;

}  // namespace physics
}  // namespace principia
