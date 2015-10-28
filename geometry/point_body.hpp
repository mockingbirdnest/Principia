#pragma once

#include <string>
#include <vector>

#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::Product;
using quantities::Quantity;
using quantities::SIUnit;

namespace geometry {

template<typename Vector>
class PointSerializer {};

template<typename Dimensions>
class PointSerializer<Quantity<Dimensions>> {
 public:
  using Vector = Quantity<Dimensions>;
  static void WriteToMessage(Vector const& coordinates,
                             not_null<serialization::Point*> const message) {
    coordinates.WriteToMessage(message->mutable_scalar());
  }

  static Vector ReadFromMessage(serialization::Point const& message) {
    CHECK(message.has_scalar());
    return Vector::ReadFromMessage(message.scalar());
  }
};

template<typename Scalar, typename Frame, int rank>
class PointSerializer<Multivector<Scalar, Frame, rank>> {
 public:
  using Vector = Multivector<Scalar, Frame, rank>;
  static void WriteToMessage(
      Vector const& coordinates,
      not_null<serialization::Point*> const message) {
    coordinates.WriteToMessage(message->mutable_multivector());
  }

  static Vector ReadFromMessage(serialization::Point const& message) {
    CHECK(message.has_multivector());
    return Vector::ReadFromMessage(message.multivector());
  }
};

template<typename Vector>
Point<Vector>::Point() {}

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
void Point<Vector>::WriteToMessage(
    not_null<serialization::Point*> const message) const {
  PointSerializer<Vector>::WriteToMessage(coordinates_, message);
}

template<typename Vector>
Point<Vector> Point<Vector>::ReadFromMessage(
    serialization::Point const& message) {
  return Point(PointSerializer<Vector>::ReadFromMessage(message));
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
std::string DebugString(Point<Vector> const& point) {
  using quantities::DebugString;
  return DebugString(point.coordinates_);
}

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Point<Vector> const& point) {
  return out << DebugString(point);
}

template<typename Vector, typename Weight>
void BarycentreCalculator<Point<Vector>, Weight>::Add(
    Point<Vector> const& point,
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

template<typename Vector, typename Weight>
Point<Vector> BarycentreCalculator<Point<Vector>, Weight>::Get() const {
  CHECK(!empty_) << "Empty BarycentreCalculator";
  return Point<Vector>(weighted_sum_ / weight_);
}

}  // namespace geometry
}  // namespace principia
