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

template<typename T1, typename T2>
class Pair {
 public:
  Pair(T1 const& t1, T2 const& t2);

  // This odd template "hides" the operator unless it is actually needed.
  // That's necessary because if T1 and T2 are vector spaces, we would end up
  // with two operator- with the same profile.
  template<typename U1 = T1, typename U2 = T2>
  typename vector_of<Pair<T1, T2>>::type operator-(
      Pair<T1, T2> const& from) const;

  Pair operator+(typename vector_of<Pair>::type const& translation) const;
  Pair operator-(typename vector_of<Pair>::type const& translation) const;

  Pair& operator+=(typename vector_of<Pair>::type const& translation);
  Pair& operator-=(typename vector_of<Pair>::type const& translation);

  bool operator==(Pair const& right) const;
  bool operator!=(Pair const& right) const;

 private:
  // This is needed so that different instantiations of Pair cannot access the
  // members.
  template<typename T1, typename T2>
  friend class Pair;

  T1 t1_;
  T2 t2_;
};

}  // namespace geometry
}  // namespace principia

#include "geometry/pair_body.hpp"
