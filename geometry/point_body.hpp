#pragma once

#include <vector>

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

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool>
operator<(Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ < right.coordinates_;
}

template<typename Vector, typename Weight>
Point<Vector> Barycentre(std::vector<Point<Vector>> const& points,
                         std::vector<Weight> const& weights) {
  CHECK_EQ(points.size(), weights.size());
  CHECK(!points.empty());
  // We need 'auto' here because we cannot easily write the type of the product
  // of a |Vector| with a |Weight|.  This is also why the loop below starts at
  // 1, as we use element 0 to get the type.
  auto weighted_sum = points[0].coordinates_ * weights[0];
  Weight weight = weights[0];
  for (size_t i = 1; i < points.size(); ++i) {
    weighted_sum += points[i].coordinates_ * weights[i];
    weight += weights[i];
  }
  return Point<Vector>(weighted_sum / weight);
}

}  // namespace geometry
}  // namespace principia
