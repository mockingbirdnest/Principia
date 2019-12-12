
#pragma once

#include "base/mappable.hpp"
#include "base/not_constructible.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/point.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace testing_utilities {
FORWARD_DECLARE_FROM(componentwise,
                     TEMPLATE(typename PairType) class,
                     ComponentwiseMatcher2Impl);
}  // namespace testing_utilities

namespace geometry {
namespace internal_pair {

using base::not_constructible;
using base::not_null;
using quantities::Product;
using quantities::Quotient;

template<typename T1, typename T2>
class Pair;

// A template to peel off the affine layer (i.e., the class Point) if any.
template<typename T>
struct vector_of : not_constructible {
  using type = T;
};

template<typename T1, typename T2>
struct vector_of<Pair<T1, T2>> : not_constructible {
  using type = Pair<typename vector_of<T1>::type, typename vector_of<T2>::type>;
};

template<typename T>
struct vector_of<Point<T>> : not_constructible {
  using type = T;
};

// A template to enable declarations on affine pairs (i.e., when one of the
// components is a Point).
template<typename T>
struct enable_if_affine : not_constructible {};

template<typename T1, typename T2>
struct enable_if_affine<Pair<Point<T1>, T2>> : not_constructible {
  using type = Pair<Point<T1>, T2>;
};

template<typename T1, typename T2>
struct enable_if_affine<Pair<T1, Point<T2>>> : not_constructible {
  using type = Pair<T1, Point<T2>>;
};

template<typename T1, typename T2>
struct enable_if_affine<Pair<Point<T1>, Point<T2>>> : not_constructible {
  using type = Pair<Point<T1>, Point<T2>>;
};

// A template to enable declarations on vector pairs (i.e., when none of the
// components is a Point).
template<typename T, typename U = T>
struct enable_if_vector : not_constructible {
  using type = U;
};

template<typename T1, typename T2>
struct enable_if_vector<Pair<Point<T1>, T2>> : not_constructible {};

template<typename T1, typename T2>
struct enable_if_vector<Pair<T1, Point<T2>>> :not_constructible {};

template<typename T1, typename T2>
struct enable_if_vector<Pair<Point<T1>, Point<T2>>> : not_constructible {};

// This class represents a pair of two values which can be members of an affine
// space (i.e., Points) or of a vector space (such as double, Quantity, Vector,
// Bivector or Trivector).  Only the operations that make sense are defined,
// depending on the nature of the parameters T1 and T2.
template<typename T1, typename T2>
class Pair {
 public:
  Pair(T1 const& t1, T2 const& t2);
  virtual ~Pair() = default;

  Pair operator+(typename vector_of<Pair>::type const& right) const;
  Pair operator-(typename vector_of<Pair>::type const& right) const;

  Pair& operator+=(typename vector_of<Pair>::type const& right);
  Pair& operator-=(typename vector_of<Pair>::type const& right);

  bool operator==(Pair const& right) const;
  bool operator!=(Pair const& right) const;

  static constexpr bool is_serializable = base::is_serializable_v<T1> &&
                                          base::is_serializable_v<T2>;

  void WriteToMessage(not_null<serialization::Pair*> message) const;
  template<typename = std::enable_if_t<is_serializable>>
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
  friend class geometry::BarycentreCalculator;

  // This is needed to make Pair mappable.
  template<typename Functor, typename T, typename>
  friend struct base::Mappable;

  // This is needed for testing.
  template<typename PairType>
  friend class testing_utilities::internal_componentwise::
      ComponentwiseMatcher2Impl;

  template<typename U1, typename U2>
  friend typename vector_of<Pair<U1, U2>>::type operator-(
      typename enable_if_affine<Pair<U1, U2>>::type const& left,
      Pair<U1, U2> const& right);

  template<typename U1, typename U2>
  friend typename enable_if_vector<Pair<U1, U2>>::type operator+(
      Pair<U1, U2> const& right);

  template<typename U1, typename U2>
  friend typename enable_if_vector<Pair<U1, U2>>::type operator-(
      Pair<U1, U2> const& right);

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
  friend typename enable_if_vector<Pair<U1, U2>>::type& operator*=(
      Pair<U1, U2>& left,
      double right);

  template<typename U1, typename U2>
  friend typename enable_if_vector<Pair<U1, U2>>::type& operator/=(
      Pair<U1, U2>& left,
      double right);

  template<typename U1, typename U2>
  friend std::ostream& operator<<(std::ostream& out, Pair<U1, U2> const& pair);
};

template<typename T1, typename T2>
typename vector_of<Pair<T1, T2>>::type operator-(
    typename enable_if_affine<Pair<T1, T2>>::type const& left,
    Pair<T1, T2> const& right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator+(
    Pair<T1, T2> const& right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator-(
    Pair<T1, T2> const& right);

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
typename enable_if_vector<Pair<T1, T2>>::type& operator*=(Pair<T1, T2>& left,
                                                          double right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type& operator/=(Pair<T1, T2>& left,
                                                          double right);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair);

}  // namespace internal_pair

using internal_pair::enable_if_vector;
using internal_pair::Pair;
using internal_pair::vector_of;

// Specialize BarycentreCalculator to make it applicable to Pairs.
namespace internal_barycentre_calculator {

template<typename T1, typename T2, typename Weight>
class BarycentreCalculator<Pair<T1, T2>, Weight> final {
 public:
  BarycentreCalculator() = default;

  void Add(Pair<T1, T2> const& pair, Weight const& weight);
  Pair<T1, T2> Get() const;

  Weight const& weight() const;

 private:
  bool empty_ = true;
  Product<typename vector_of<T1>::type, Weight> t1_weighted_sum_;
  Product<typename vector_of<T2>::type, Weight> t2_weighted_sum_;
  Weight weight_;

  // We need reference values to convert points into vectors, if needed.  We
  // pick default-constructed objects as they don't introduce any inaccuracies
  // in the computations.
  static T1 const reference_t1_;
  static T2 const reference_t2_;
};

}  // namespace internal_barycentre_calculator
}  // namespace geometry

// Reopen the base namespace to make Pairs of vectors mappable.
namespace base {

template<typename Functor, typename T1, typename T2>
class Mappable<Functor,
               geometry::Pair<T1, T2>,
               typename geometry::enable_if_vector<
                   geometry::Pair<T1, T2>, void>::type> {
 public:
  using type = geometry::Pair<
                   decltype(std::declval<Functor>()(std::declval<T1>())),
                   decltype(std::declval<Functor>()(std::declval<T2>()))>;

  static type Do(Functor const& functor,
                 geometry::Pair<T1, T2> const& pair);
};

}  // namespace base
}  // namespace principia

#include "geometry/pair_body.hpp"
