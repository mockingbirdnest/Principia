#pragma once

#include <string>

#include "base/concepts.hpp"
#include "base/mappable.hpp"  // 🧙 For base::_mappable::internal.
#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/barycentre_calculator.hpp"  // 🧙 For friendship.
#include "geometry/space.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE(TEMPLATE(typename PairType) class,
                ComponentwiseMatcher2Impl,
                FROM(componentwise));
}  // namespace testing_utilities

namespace geometry {
namespace _pair {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_space;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;

template<typename T1, typename T2>
class Pair;

// A template to peel off the affine layer, if any, in a Pair.
template<typename T>
struct vector_of;

template<typename T1, typename T2>
struct vector_of<Pair<T1, T2>> : not_constructible {
  using type = Pair<Difference<T1>, Difference<T2>>;
};

template<typename T>
using vector_of_t = typename vector_of<T>::type;

// A template to enable declarations on affine pairs (i.e., when one of the
// components is not a vector).
template<typename T, typename U = T, typename = void>
struct enable_if_affine;

template<typename T1, typename T2, typename U>
struct enable_if_affine<
    Pair<T1, T2>, U,
    std::enable_if_t<!(additive_group<T1> && additive_group<T2>)>>
    : not_constructible {
  using type = U;
};

template<typename T>
using enable_if_affine_t = typename enable_if_affine<T>::type;

template<typename T>
struct is_position : std::false_type {};
template<typename Frame>
struct is_position<Position<Frame>> : std::true_type {};
template<typename T>
constexpr bool is_position_v = is_position<T>::value;

template<typename T>
struct is_displacement : std::false_type {};
template<typename Frame>
struct is_displacement<Displacement<Frame>> : std::true_type {};
template<typename T>
constexpr bool is_displacement_v = is_displacement<T>::value;

template<typename T>
struct is_velocity : std::false_type {};
template<typename Frame>
struct is_velocity<Velocity<Frame>> : std::true_type {};
template<typename T>
constexpr bool is_velocity_v = is_velocity<T>::value;

// A template to enable declarations on vector pairs (i.e., when both of the
// components are vectors).
template<typename T, typename U = T, typename = void>
struct enable_if_vector;

template<typename T1, typename T2, typename U>
struct enable_if_vector<
    Pair<T1, T2>, U,
    std::enable_if_t<additive_group<T1> && additive_group<T2>>>
    : not_constructible {
  using type = U;
};

template<typename T, typename U = T>
using enable_if_vector_t = typename enable_if_vector<T, U>::type;

// This class represents a pair of two values which can be members of an affine
// space (i.e., Points) or of a vector space (such as double, Quantity, Vector,
// Bivector or Trivector).  Only the operations that make sense are defined,
// depending on the nature of the parameters T1 and T2.
template<typename T1, typename T2>
class Pair final {
 public:
  Pair()
    requires std::default_initializable<T1> && std::default_initializable<T2> =
      default;
  Pair(T1 const& t1, T2 const& t2);

  friend bool operator==(Pair const& left, Pair const& right) = default;
  friend bool operator!=(Pair const& left, Pair const& right) = default;

  Pair operator+(vector_of_t<Pair> const& right) const;
  Pair operator-(vector_of_t<Pair> const& right) const;

  Pair& operator+=(vector_of_t<Pair> const& right);
  Pair& operator-=(vector_of_t<Pair> const& right);
  template<typename U1 = T1, typename U2 = T2>
  enable_if_vector_t<Pair<U1, U2>>& operator*=(double right);
  template<typename U1 = T1, typename U2 = T2>
  enable_if_vector_t<Pair<U1, U2>>& operator/=(double right);

  T1 const& position() const
    requires is_position_v<T1>;
  T1 const& displacement() const
    requires is_displacement_v<T1>;

  T2 const& velocity() const
    requires is_velocity_v<T2>;

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  static Pair ReadFromMessage(serialization::Pair const& message)
    requires serializable<T1> && serializable<T2>;

 private:
  T1 t1_;
  T2 t2_;

  // This is needed so that different instantiations of Pair can access the
  // members.
  template<typename U1, typename U2>
  friend class Pair;

  // This is needed to make Pair mappable.
  template<typename Functor, typename T, typename>
  friend struct base::_mappable::internal::Mappable;

  // This is needed to implement matchers for `Pair`.
  template<typename PairType>
  friend class testing_utilities::_componentwise::internal::
      ComponentwiseMatcher2Impl;

  template<typename U1, typename U2>
  friend vector_of_t<Pair<U1, U2>> operator-(
      enable_if_affine_t<Pair<U1, U2>> const& left,
      Pair<U1, U2> const& right);

  template<typename U1, typename U2>
  friend enable_if_vector_t<Pair<U1, U2>> operator+(Pair<U1, U2> const& right);

  template<typename U1, typename U2>
  friend enable_if_vector_t<Pair<U1, U2>> operator-(Pair<U1, U2> const& right);

  template<typename Scalar, typename U1, typename U2>
  friend typename enable_if_vector<
      Pair<U1, U2>,
      Pair<Product<Scalar, U1>, Product<Scalar, U2>>>::type
  operator*(Scalar left, Pair<U1, U2> const& right);

  template<typename Scalar, typename U1, typename U2>
  friend typename enable_if_vector<
      Pair<U1, U2>,
      Pair<Product<U1, Scalar>, Product<U2, Scalar>>>::type
  operator*(Pair<U1, U2> const& left, Scalar right);

  template<typename Scalar, typename U1, typename U2>
  friend typename enable_if_vector<
      Pair<U1, U2>,
      Pair<Quotient<U1, Scalar>, Quotient<U2, Scalar>>>::type
  operator/(Pair<U1, U2> const& left, Scalar right);

  template<typename U1, typename U2>
  friend std::string DebugString(Pair<U1, U2> const& pair);

  template<typename U1, typename U2>
  friend std::ostream& operator<<(std::ostream& out, Pair<U1, U2> const& pair);
};

template<typename T1, typename T2>
vector_of_t<Pair<T1, T2>> operator-(
    enable_if_affine_t<Pair<T1, T2>> const& left,
    Pair<T1, T2> const& right);

template<typename T1, typename T2>
enable_if_vector_t<Pair<T1, T2>> operator+(Pair<T1, T2> const& right);

template<typename T1, typename T2>
enable_if_vector_t<Pair<T1, T2>> operator-(Pair<T1, T2> const& right);

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Product<Scalar, T1>, Product<Scalar, T2>>>::type
operator*(Scalar left, Pair<T1, T2> const& right);

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Product<T1, Scalar>, Product<T2, Scalar>>>::type
operator*(Pair<T1, T2> const& left, Scalar right);

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<Quotient<T1, Scalar>, Quotient<T2, Scalar>>>::type
operator/(Pair<T1, T2> const& left, Scalar right);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair);

}  // namespace internal

using internal::enable_if_vector_t;
using internal::Pair;
using internal::vector_of;

}  // namespace _pair
}  // namespace geometry

// Reopen the base namespace to make Pairs of vectors mappable.
namespace base {
namespace _mappable {
namespace internal {

using namespace principia::geometry::_pair;
using namespace principia::quantities::_named_quantities;

template<typename Functor, typename T1, typename T2>
class Mappable<Functor,
               Pair<T1, T2>,
               enable_if_vector_t<Pair<T1, T2>, void>> {
 public:
  using type = Pair<decltype(std::declval<Functor>()(std::declval<T1>())),
                    decltype(std::declval<Functor>()(std::declval<T2>()))>;

  static type Do(Functor const& functor, Pair<T1, T2> const& pair);
};

}  // namespace internal
}  // namespace _mappable
}  // namespace base
}  // namespace principia

#include "geometry/pair_body.hpp"
