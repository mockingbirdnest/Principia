#pragma once

#include "base/algebra.hpp"
#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _space {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Frame>
using Displacement = Vector<Length, Frame>;

template<typename Frame>
using Position = Point<Displacement<Frame>>;

template<typename Frame>
using Velocity = Vector<Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<AngularFrequency, Frame>;

}  // namespace internal

using internal::AngularVelocity;
using internal::Displacement;
using internal::Position;
using internal::Velocity;

}  // namespace _space

namespace _direct_sum {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_not_null;
using namespace principia::base::_tags;
using namespace principia::geometry::_space;

// Specializations for use by `physics`.  They are declared here to make sure
// that code that uses `_space` sees the specializations.

template<typename Frame>
class DirectSum<Position<Frame>, Velocity<Frame>> {
 public:
  constexpr DirectSum() = default;
  constexpr explicit DirectSum(uninitialized_t);

  constexpr DirectSum(Position<Frame> const& position,
                      Velocity<Frame> const& velocity);

  bool operator==(DirectSum const&) const = default;

  DirectSum& operator+=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);
  DirectSum& operator-=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  static DirectSum<Position<Frame>, Velocity<Frame>> ReadFromMessage(
      serialization::Pair const& message);

 private:
  Position<Frame> position_;
  Velocity<Frame> velocity_;
};

}  // namespace internal
}  // namespace _direct_sum
}  // namespace geometry
}  // namespace principia
