#pragma once

#include "geometry/pair.hpp"

namespace principia {
namespace geometry {

template<typename T1, typename T2>
Pair<T1, T2>::Pair(T1 const& t1, T2 const& t2)
    : t1_(t1),
      t2_(t2) {}

template<typename T1, typename T2>
template<typename U1, typename U2>
typename vector_of<Pair<T1, T2>>::type Pair<T1, T2>::operator-(
    Pair<T1, T2> const& from) const {
  return typename vector_of<Pair<T1, T2>>::type(t1_ - from.t1_, t2_ - from.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator+(
    typename vector_of<Pair>::type const& translation) const {
  return Pair<T1, T2>(t1_ + translation.t1_, t2_ + translation.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2> Pair<T1, T2>::operator-(
    typename vector_of<Pair>::type const& translation) const {
  return Pair<T1, T2>(t1_ - translation.t1_, t2_ - translation.t2_);
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator+=(
    typename vector_of<Pair>::type const& translation) {
  t1_ += translation.t1_;
  t2_ += translation.t2_;
  return *this;
}

template<typename T1, typename T2>
Pair<T1, T2>& Pair<T1, T2>::operator-=(
    typename vector_of<Pair>::type const& translation) {
  t1_ -= translation.t1_;
  t2_ -= translation.t2_;
  return *this;
}

template<typename T1, typename T2>
bool Pair<T1, T2>::operator==(Pair const& right) const {
  return t1_ == right.t1_ && t2_ == right.t2_;
}

template<typename T1, typename T2>
bool Pair<T1, T2>::operator!=(Pair const& right) const {
  return t1_ != right.t1_ || t2_ != right.t2_;
}

}  // namespace geometry
}  // namespace principia
