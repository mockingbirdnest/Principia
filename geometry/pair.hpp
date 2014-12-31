#pragma once

#include "geometry/point.hpp"

namespace principia {
namespace geometry {

template<typename T1, typename T2>
class Pair;

//TODO(phl): Style for these auxiliary templates?

template<typename T>
class VectorOf {
 public:
  using type = T;
};

template<typename T1, typename T2>
class VectorOf<Pair<T1, T2>> {
 public:
  using type = Pair<typename VectorOf<T1>::type, typename VectorOf<T2>::type>;
};

template<typename T>
class VectorOf<Point<T>> {
 public:
  using type = T;
};

template<typename T1, typename T2>
class Pair {
 public:
  Pair(T1 const& t1, T2 const& t2);

  typename VectorOf<Pair>::type operator-(Pair const& from) const;

  Pair operator+(typename VectorOf<Pair>::type const& translation) const;
  Pair operator-(typename VectorOf<Pair>::type const& translation) const;

  Pair& operator+=(typename VectorOf<Pair>::type const& translation);
  Pair& operator-=(typename VectorOf<Pair>::type const& translation);

  bool operator==(Pair const& right) const;
  bool operator!=(Pair const& right) const;

 private:
  T1 t1_;
  T2 t2_;
};

}  // namespace geometry
}  // namespace principia

#include "geometry/pair_body.hpp"
