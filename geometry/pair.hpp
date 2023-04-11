#pragma once

#include "base/mappable.hpp"
#include "base/not_constructible.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/point.hpp"
#include "geometry/traits.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE_FROM(componentwise,
                     TEMPLATE(typename PairType) class,
                     ComponentwiseMatcher2Impl);
}  // namespace testing_utilities

namespace geometry {
namespace _pair {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_traits;
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
    std::enable_if_t<!std::conjunction_v<is_vector<T1>, is_vector<T2>>>>
    : not_constructible {
  using type = U;
};

template<typename T>
using enable_if_affine_t = typename enable_if_affine<T>::type;

// A template to enable declarations on vector pairs (i.e., when both of the
// components are vectors).
template<typename T, typename U = T, typename = void>
struct enable_if_vector;

template<typename T1, typename T2, typename U>
struct enable_if_vector<
    Pair<T1, T2>, U,
    std::enable_if_t<std::conjunction_v<is_vector<T1>, is_vector<T2>>>>
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
class Pair {
 public:
  Pair(T1 const& t1, T2 const& t2);
  virtual ~Pair() = default;

  Pair operator+(vector_of_t<Pair> const& right) const;
  Pair operator-(vector_of_t<Pair> const& right) const;

  Pair& operator+=(vector_of_t<Pair> const& right);
  Pair& operator-=(vector_of_t<Pair> const& right);

  bool operator==(Pair const& right) const;
  bool operator!=(Pair const& right) const;

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  template<typename U1 = T1,
           typename U2 = T2,
           typename = std::enable_if_t<is_serializable_v<U1> &&
                                       is_serializable_v<U2>>>
  static Pair ReadFromMessage(serialization::Pair const& message);

 protected:
  // The subclasses can access the members directly to implement accessors.
  T1 t1_;
  T2 t2_;

 private:
  // This is needed so that different instantiations of Pair can access the
  // members.
  template<typename U1, typename U2>
  friend class Pair;

  // This is needed to specialize BarycentreCalculator.
  template<typename V, typename S>
  friend class _barycentre_calculator::BarycentreCalculator;

  // This is needed to make Pair mappable.
  template<typename Functor, typename T, typename>
  friend struct base::_mappable::internal::Mappable;

  // This is needed for testing.
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
  friend enable_if_vector_t<Pair<U1, U2>>& operator*=(Pair<U1, U2>& left,
                                                      double right);

  template<typename U1, typename U2>
  friend enable_if_vector_t<Pair<U1, U2>>& operator/=(Pair<U1, U2>& left,
                                                      double right);

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
enable_if_vector_t<Pair<T1, T2>>& operator*=(Pair<T1, T2>& left, double right);

template<typename T1, typename T2>
enable_if_vector_t<Pair<T1, T2>>& operator/=(Pair<T1, T2>& left, double right);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair);

}  // namespace internal

using internal::enable_if_vector_t;
using internal::Pair;
using internal::vector_of;

}  // namespace _pair

// Specialize BarycentreCalculator to make it applicable to Pairs.
namespace _barycentre_calculator {
namespace internal {

using namespace principia::geometry::_pair;

template<typename T1, typename T2, typename Weight>
class BarycentreCalculator<Pair<T1, T2>, Weight> final {
 public:
  BarycentreCalculator() = default;

  void Add(Pair<T1, T2> const& pair, Weight const& weight);
  Pair<T1, T2> Get() const;

  Weight const& weight() const;

 private:
  bool empty_ = true;
  Product<Difference<T1>, Weight> t1_weighted_sum_;
  Product<Difference<T2>, Weight> t2_weighted_sum_;
  Weight weight_;

  // We need reference values to convert points into vectors, if needed.  We
  // pick default-constructed objects as they don't introduce any inaccuracies
  // in the computations.
  static T1 const reference_t1_;
  static T2 const reference_t2_;
};

}  // namespace internal
}  // namespace _barycentre_calculator
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
