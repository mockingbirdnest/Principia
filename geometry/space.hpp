#pragma once

#include <cstdint>

#include "base/algebra.hpp"
#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "geometry/direct_sum.hpp"
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
using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_space;

// Specializations for use by `physics`.  They are declared in this file to make
// sure that code that uses `_space` sees the specializations.

template<typename Frame>
class DirectSum<Position<Frame>, Velocity<Frame>> {
 public:
  constexpr DirectSum() = default;

  // This constructor actually initializes the members.
  constexpr explicit DirectSum(uninitialized_t);

  constexpr DirectSum(Position<Frame> const& position,
                      Velocity<Frame> const& velocity);

  // Visible by ADL.
  template<std::size_t i>
  friend constexpr auto const& get(DirectSum const& direct_sum) {
    if constexpr (i == 0) {
      return direct_sum.position();
    } else if constexpr (i == 1) {
      return direct_sum.velocity();
    } else {
      static_assert(
          i < 2, "Index out of bounds in get<DirectSum<Position, Velocity>>");
    }
  }

  //TODO(phl)Can we hide this?
  template<std::size_t i>
  friend constexpr auto& get(DirectSum& direct_sum) {
    if constexpr (i == 0) {
      return direct_sum.position();
    } else if constexpr (i == 1) {
      return direct_sum.velocity();
    } else {
      static_assert(
          i < 2, "Index out of bounds in get<DirectSum<Position, Velocity>>");
    }
  }

  bool operator==(DirectSum const&) const = default;

  DirectSum& operator+=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);
  DirectSum& operator-=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);

  template<typename Self>
  constexpr auto&& position(this Self&& self);
  template<typename Self>
  constexpr auto&& velocity(this Self&& self);

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  static DirectSum ReadFromMessage(serialization::Pair const& message);

 private:
  Position<Frame> position_;
  Velocity<Frame> velocity_;
};

template<typename Frame>
class DirectSum<Displacement<Frame>, Velocity<Frame>> {
 public:
  constexpr DirectSum() = default;

  // This constructor actually initializes the members.
  constexpr explicit DirectSum(uninitialized_t);

  constexpr DirectSum(Displacement<Frame> const& displacement,
                      Velocity<Frame> const& velocity);

  // Visible by ADL.
  template<std::size_t i>
  friend constexpr auto const& get(DirectSum const& direct_sum) {
    if constexpr (i == 0) {
      return direct_sum.displacement();
    } else if constexpr (i == 1) {
      return direct_sum.velocity();
    } else {
      static_assert(
          i < 2,
          "Index out of bounds in get<DirectSum<Displacement, Velocity>>");
    }
  }

  //TODO(phl)Can we hide this?
  template<std::size_t i>
  friend constexpr auto& get(DirectSum& direct_sum) {
    if constexpr (i == 0) {
      return direct_sum.displacement();
    } else if constexpr (i == 1) {
      return direct_sum.velocity();
    } else {
      static_assert(
          i < 2,
          "Index out of bounds in get<DirectSum<Displacement, Velocity>>");
    }
  }

  bool operator==(DirectSum const&) const = default;

  DirectSum& operator+=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);
  DirectSum& operator-=(
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& right);

  template<ring Scalar>
  DirectSum& operator*=(Scalar const& right)
    requires(module<Displacement<Frame>, Scalar> &&
             module<Velocity<Frame>, Scalar>);
  template<field Scalar>
  DirectSum& operator/=(Scalar const& right)
    requires(vector_space<Displacement<Frame>, Scalar> &&
             vector_space<Velocity<Frame>, Scalar>);

  template<typename Self>
  constexpr auto&& displacement(this Self&& self);
  template<typename Self>
  constexpr auto&& velocity(this Self&& self);

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  static DirectSum ReadFromMessage(serialization::Pair const& message);

 private:
  Displacement<Frame> displacement_;
  Velocity<Frame> velocity_;
};

template<std::size_t i, typename Frame>
constexpr auto const& get(
    DirectSum<Position<Frame>, Velocity<Frame>> const& direct_sum);

template<std::size_t i, typename Frame>
constexpr auto const& get(
    DirectSum<Displacement<Frame>, Velocity<Frame>> const& direct_sum);

}  // namespace internal
}  // namespace _direct_sum
}  // namespace geometry

// Reopen the base namespace to make DirectSum<Displacement, Velocity> mappable.
namespace base {
namespace _mappable {
namespace internal {

using namespace principia::geometry::_direct_sum;
using namespace principia::geometry::_space;

template<typename Functor, typename Frame>
class Mappable<Functor, DirectSum<Displacement<Frame>, Velocity<Frame>>> {
 public:
  using type = DirectSum<
      decltype(std::declval<Functor>()(std::declval<Displacement<Frame>>())),
      decltype(std::declval<Functor>()(std::declval<Velocity<Frame>>()))>;

  static type Do(
      Functor const& functor,
      DirectSum<Displacement<Frame>, Velocity<Frame>> const& direct_sum);
};

}  // namespace internal
}  // namespace _mappable
}  // namespace base

}  // namespace principia

#include "geometry/space_body.hpp"
