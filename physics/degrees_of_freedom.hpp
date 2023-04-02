#pragma once

#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/pair.hpp"
#include "geometry/point.hpp"
#include "geometry/space.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _degrees_of_freedom {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_traits;
using namespace principia::geometry::_pair;
using namespace principia::geometry::_space;

// This class is analogous to the pair which is its base class, except that it
// exports properly-named selectors.  It is implicitly convertible in both
// directions, so clients can generally ignore the difference.  Note however
// that creating a DegreesOfFreedom involves a copy so clients might want to use
// the base type (probably declared as |auto|) when they don't need to access
// the members.
template<typename Frame>
class DegreesOfFreedom : public Pair<Position<Frame>, Velocity<Frame>> {
 public:
  DegreesOfFreedom(Position<Frame> const& position,
                   Velocity<Frame> const& velocity);

  // Not explicit, the point of this constructor is to convert implicitly.
  DegreesOfFreedom(
      Pair<Position<Frame>,
           Velocity<Frame>> const& base);  // NOLINT(runtime/explicit)

  template<typename F = Frame,
           typename = std::enable_if_t<is_serializable_v<F>>>
  static DegreesOfFreedom ReadFromMessage(serialization::Pair const& message);

  Position<Frame> const& position() const;
  Velocity<Frame> const& velocity() const;
};

// This class is analogous to the vector class underlying DegreesOfFreedom,
// except that it exports properly-named selectors.  The same comments as above
// apply.
template<typename Frame>
class RelativeDegreesOfFreedom
    : public Pair<Displacement<Frame>, Velocity<Frame>> {
 public:
  RelativeDegreesOfFreedom() = default;

  RelativeDegreesOfFreedom(Displacement<Frame> const& displacement,
                           Velocity<Frame> const& velocity);

  // Not explicit, the point of this constructor is to convert implicitly.
  RelativeDegreesOfFreedom(
      Pair<Displacement<Frame>,
           Velocity<Frame>> const& base);  // NOLINT(runtime/explicit)

  Displacement<Frame> const& displacement() const;
  Velocity<Frame> const& velocity() const;
};

template<typename Frame>
std::string DebugString(DegreesOfFreedom<Frame> const& degrees_of_freedom);

template<typename Frame>
std::string DebugString(
    RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom);

template<typename Frame>
std::ostream& operator<<(std::ostream& out,
                         DegreesOfFreedom<Frame> const& degrees_of_freedom);

template<typename Frame>
std::ostream& operator<<(
    std::ostream& out,
    RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom);

}  // namespace internal

using internal::DegreesOfFreedom;
using internal::RelativeDegreesOfFreedom;

}  // namespace _degrees_of_freedom
}  // namespace physics

// Reopen the base namespace to make RelativeDegreesOfFreedom mappable.
namespace base {
namespace _mappable {
namespace internal {

using namespace principia::geometry::_pair;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;

template<typename Functor, typename Frame>
struct Mappable<Functor, RelativeDegreesOfFreedom<Frame>>
    : not_constructible {
  using type = Pair<decltype(std::declval<Functor>()(
                        std::declval<Displacement<Frame>>())),
                    decltype(std::declval<Functor>()(
                        std::declval<Velocity<Frame>>()))>;

  static type Do(Functor const& functor,
                 RelativeDegreesOfFreedom<Frame> const& relative);
};

}  // namespace internal
}  // namespace _mappable
}  // namespace base

// Reopen the geometry namespace to make BarycentreCalculator applicable to
// degrees of freedom.
namespace geometry {
namespace _barycentre_calculator {
namespace internal {

using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;

template<typename Frame, typename Weight>
class BarycentreCalculator<DegreesOfFreedom<Frame>, Weight> final {
 public:
  BarycentreCalculator() = default;

  void Add(DegreesOfFreedom<Frame> const& degrees_of_freedom,
           Weight const& weight);
  DegreesOfFreedom<Frame> Get() const;

  Weight const& weight() const;

 private:
  BarycentreCalculator<Pair<Position<Frame>, Velocity<Frame>>, Weight>
      implementation_;
};

template<typename Frame, typename Weight>
class BarycentreCalculator<RelativeDegreesOfFreedom<Frame>, Weight>
    final {
 public:
  BarycentreCalculator() = default;

  void Add(RelativeDegreesOfFreedom<Frame> const& relative_degrees_of_freedom,
           Weight const& weight);
  RelativeDegreesOfFreedom<Frame> Get() const;

  Weight const& weight() const;

 private:
  BarycentreCalculator<Pair<Displacement<Frame>, Velocity<Frame>>, Weight>
      implementation_;
};

}  // namespace internal
}  // namespace _barycentre_calculator
}  // namespace geometry
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_degrees_of_freedom;
}  // namespace principia::physics

#include "physics/degrees_of_freedom_body.hpp"
