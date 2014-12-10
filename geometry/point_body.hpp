#pragma once

#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

using principia::quantities::Product;
using principia::quantities::SIUnit;

namespace principia {
namespace geometry {

template<typename Vector>
Point<Vector>::Point() {}

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
template<typename Weight>
void Point<Vector>::BarycentreCalculator<Weight>::Add(Point const& point,
                                                      Weight const& weight) {
  if (empty_) {
    weighted_sum_ = point.coordinates_ * weight;
    weight_ = weight;
    empty_ = false;
  } else {
    weighted_sum_ += point.coordinates_ * weight;
    weight_ += weight;
  }
}

template<typename Vector>
template<typename Weight>
Point<Vector> const Point<Vector>::BarycentreCalculator<Weight>::Get() const {
  CHECK(!empty_) << "Empty BarycentreCalculator";
  return Point<Vector>(weighted_sum_ / weight_);
}

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point) {
  return Point<Vector>(translation + point.coordinates_);
}

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator<(
    Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ < right.coordinates_;
}

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator<=(
    Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ <= right.coordinates_;
}

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator>=(
    Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ >= right.coordinates_;
}

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator>(
    Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ > right.coordinates_;
}

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Point<Vector> const& point) {
  return out << point.coordinates_;
}

template<typename Vector, typename Weight>
Point<Vector> Barycentre(std::vector<Point<Vector>> const& points,
                         std::vector<Weight> const& weights) {
  CHECK_EQ(points.size(), weights.size());
  CHECK(!points.empty());
  Point<Vector>::BarycentreCalculator<Weight> calculator;
  for (size_t i = 1; i < points.size(); ++i) {
    calculator.Add(points[i], weights[i]);
  }
  return calculator.Get();
}

}  // namespace geometry
}  // namespace principia
