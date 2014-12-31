#pragma once

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

// A template to enable declarations on vector pairs (i.e., none of the
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

 private:
  T1 t1_;
  T2 t2_;

  // This is needed so that different instantiations of Pair cannot access the
  // members.
  template<typename T1, typename T2>
  friend class Pair;

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

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator*(
    Pair<T1, T2> const& left,
    double const right);

template<typename T1, typename T2>
typename enable_if_vector<Pair<T1, T2>>::type operator/(
    Pair<T1, T2> const& left,
    double const right);

}  // namespace geometry
}  // namespace principia

#include "geometry/pair_body.hpp"
