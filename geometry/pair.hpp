#pragma once

#include "base/mappable.hpp"
#include "geometry/point.hpp"

namespace principia {
namespace geometry {

template<typename T1, typename T2>
class Pair;

// A template to peel off the affine layer (i.e., the class Point) if any.
template<typename T>
class vector_of {
 public:
  using type = T;
};

template<typename T1, typename T2>
class vector_of<Pair<T1, T2>> {
 public:
  using type = Pair<typename vector_of<T1>::type, typename vector_of<T2>::type>;
};

template<typename T>
class vector_of<Point<T>> {
 public:
  using type = T;
};

// A template to enable declarations on affine pairs (i.e., when one of the
// components is a Point).
template<typename T>
class enable_if_affine {};

template<typename T1, typename T2>
class enable_if_affine<Pair<Point<T1>, T2>> {
 public:
  using type = Pair<Point<T1>, T2>;
};

template<typename T1, typename T2>
class enable_if_affine<Pair<T1, Point<T2>>> {
 public:
  using type = Pair<T1, Point<T2>>;
};

template<typename T1, typename T2>
class enable_if_affine<Pair<Point<T1>, Point<T2>>> {
 public:
  using type = Pair<Point<T1>, Point<T2>>;
};

// A template to enable declarations on vector pairs (i.e., when none of the
// components is a Point).
template<typename T, typename U = T>
class enable_if_vector {
 public:
  using type = U;
};

template<typename T1, typename T2>
class enable_if_vector<Pair<Point<T1>, T2>> {};

template<typename T1, typename T2>
class enable_if_vector<Pair<T1, Point<T2>>> {};

template<typename T1, typename T2>
class enable_if_vector<Pair<Point<T1>, Point<T2>>> {};

// This class represents a pair of two values which can be members of an affine
// space (i.e., Points) or of a vector space (such as double, Quantity, Vector,
// Bivector or Trivector).  Only the operations that make sense are defined,
// depending on the nature of the parameters T1 and T2.
template<typename T1, typename T2>
class Pair {
 public:
  Pair(T1 const& t1, T2 const& t2);

  Pair operator+(typename vector_of<Pair>::type const& right) const;
  Pair operator-(typename vector_of<Pair>::type const& right) const;

  Pair& operator+=(typename vector_of<Pair>::type const& right);
  Pair& operator-=(typename vector_of<Pair>::type const& right);

  bool operator==(Pair const& right) const;
  bool operator!=(Pair const& right) const;

  template<typename Weight>
  class BarycentreCalculator {
   public:
    BarycentreCalculator() = default;
    ~BarycentreCalculator() = default;

    void Add(Pair const& pair, Weight const& weight);
    Pair const Get() const;

   private:
    bool empty_ = true;
    decltype(std::declval<typename vector_of<T1>::type>() *
             std::declval<Weight>()) t1_weighted_sum_;
    decltype(std::declval<typename vector_of<T2>::type>() *
             std::declval<Weight>()) t2_weighted_sum_;
    Weight weight_;

    // We need reference values to convert points into vectors, if needed.  We
    // pick default-constructed objects as they don't introduce any inaccuracies
    // in the computations.
    static T1 const reference_t1_;
    static T2 const reference_t2_;
  };

 protected:
  // The subclasses can access the members directly to implement accessors.
  T1 t1_;
  T2 t2_;

 private:
  // This is needed so that different instantiations of Pair can access the
  // members.
  template<typename T1, typename T2>
  friend class Pair;

  template<typename Functor, typename T, typename>
  friend class base::Mappable;

  template<typename T1, typename T2>
  friend typename vector_of<Pair<T1, T2>>::type operator-(
      typename enable_if_affine<Pair<T1, T2>>::type const& left,
      Pair<T1, T2> const& right);

  template<typename T1, typename T2>
  friend typename enable_if_vector<Pair<T1, T2>>::type operator+(
      Pair<T1, T2> const& right);

  template<typename T1, typename T2>
  friend typename enable_if_vector<Pair<T1, T2>>::type operator-(
      Pair<T1, T2> const& right);

  template<typename Scalar, typename T1, typename T2>
  friend typename enable_if_vector<
      Pair<T1, T2>,
      Pair<decltype(std::declval<Scalar>() * std::declval<T1>()),
           decltype(std::declval<Scalar>() * std::declval<T2>())>>::type
  operator*(Scalar const left, Pair<T1, T2> const& right);

  template<typename Scalar, typename T1, typename T2>
  friend typename enable_if_vector<
      Pair<T1, T2>,
      Pair<decltype(std::declval<T1>() * std::declval<Scalar>()),
           decltype(std::declval<T2>() * std::declval<Scalar>())>>::type
  operator*(Pair<T1, T2> const& left, Scalar const right);

  template<typename Scalar, typename T1, typename T2>
  friend typename enable_if_vector<
      Pair<T1, T2>,
      Pair<decltype(std::declval<T1>() * std::declval<Scalar>()),
           decltype(std::declval<T2>() * std::declval<Scalar>())>>::type
  operator/(Pair<T1, T2> const& left, Scalar const right);

  template<typename T1, typename T2>
  friend typename enable_if_vector<Pair<T1, T2>>::type& operator*=(
      Pair<T1, T2>& left,  // NOLINT(runtime/references)
      double const right);

  template<typename T1, typename T2>
  friend typename enable_if_vector<Pair<T1, T2>>::type& operator/=(
      Pair<T1, T2>& left,  // NOLINT(runtime/references)
      double const right);

  template<typename T1, typename T2>
  friend std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair);
};

// NOTE(phl): Would like to put the enable_if_affine<> on the return type, but
// this confuses MSVC.
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
    Pair<decltype(std::declval<Scalar>() * std::declval<T1>()),
         decltype(std::declval<Scalar>() * std::declval<T2>())>>::type
operator*(Scalar const left, Pair<T1, T2> const& right);

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<decltype(std::declval<T1>() * std::declval<Scalar>()),
         decltype(std::declval<T2>() * std::declval<Scalar>())>>::type
operator*(Pair<T1, T2> const& left, Scalar const right);

template<typename Scalar, typename T1, typename T2>
typename enable_if_vector<
    Pair<T1, T2>,
    Pair<decltype(std::declval<T1>() / std::declval<Scalar>()),
         decltype(std::declval<T2>() / std::declval<Scalar>())>>::type
operator/(Pair<T1, T2> const& left, Scalar const right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type& operator*=(
    Pair<T1, T2>& left,  // NOLINT(runtime/references)
    double const right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type& operator/=(
    Pair<T1, T2>& left,  // NOLINT(runtime/references)
    double const right);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, Pair<T1, T2> const& pair);

}  // namespace geometry

// Reopen the base namespace to make Pair mappable.
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
