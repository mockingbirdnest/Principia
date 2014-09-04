#pragma once

#include "glog/logging.h"
#include "quantities/quantities.hpp"

using principia::quantities::Product;
using principia::quantities::SIUnit;

namespace principia {
namespace geometry {

template<typename Vector>
Point<Vector>::Point(Vector const& coordinates) : coordinates_(coordinates) {}

template<typename Vector>
Vector Point<Vector>::operator-(Point const& from) const {
  return coordinates_ - from.coordinates_;
}

template<typename Vector>
Point<Vector> Point<Vector>::operator+(Vector const& translation) const {
  return Point<Vector>(coordinates_ + translation);
}

template<typename Vector>
Point<Vector> Point<Vector>::operator-(Vector const& translation) const {
  return Point<Vector>(coordinates_ - translation);
}

template<typename Vector>
Point<Vector>& Point<Vector>::operator+=(Vector const& translation) {
  coordinates_ += translation;
  return *this;
}

template<typename Vector>
Point<Vector>& Point<Vector>::operator-=(Vector const& translation) {
  coordinates_ -= translation;
  return *this;
}

template<typename Vector>
bool Point<Vector>::operator==(Point<Vector> const& right) const {
  return coordinates_ == right.coordinates_;
}

template<typename Vector>
bool Point<Vector>::operator!=(Point<Vector> const& right) const {
  return coordinates_ != right.coordinates_;
}

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point) {
  return Point<Vector>(translation + point.coordinates_);
}

template<typename Vector, typename Weight>
Point<Vector> Barycentre(std::vector<Point<Vector>> const& points,
                         std::vector<Weight> const& weights) {
  CHECK_EQ(points.size(), weights.size());
  Product<Vector, Weight> coordinates;
  Weight weight = 0 * SIUnit<Weight>();
  for (size_t i = 0; i < points.size(); ++i) {
    coordinates += points[i].coordinates_ * weights[i];
    weight += weights[i];
  }
  return Point<Vector>(coordinates / weight);
}

}  // namespace geometry
}  // namespace principia
