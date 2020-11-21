
#pragma once

#include "geometry/point.hpp"

#include <string>
#include <vector>

#include "base/not_constructible.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_point {

using base::not_constructible;
using quantities::Product;
using quantities::Quantity;

template<typename Vector>
struct PointSerializer : not_constructible {};

template<typename Dimensions>
struct PointSerializer<Quantity<Dimensions>> : not_constructible {
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
struct PointSerializer<Multivector<Scalar, Frame, rank>> : not_constructible {
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
constexpr Point<Vector>::Point() : coordinates_() {}

#if PRINCIPIA_COMPILER_MSVC && !__INTELLISENSE__
template<typename Vector>
constexpr Point<Vector>::Point(Point const& other)
    : coordinates_(other.coordinates_) {}

template<typename Vector>
constexpr Point<Vector>::Point(Point&& other)
    : coordinates_(std::move(other.coordinates_)) {}
#endif

template<typename Vector>
constexpr Vector Point<Vector>::operator-(Point const& from) const {
  return coordinates_ - from.coordinates_;
}

template<typename Vector>
constexpr Point<Vector> Point<Vector>::operator+(
    Vector const& translation) const {
  return Point(coordinates_ + translation);
}

template<typename Vector>
constexpr Point<Vector> Point<Vector>::operator-(
    Vector const& translation) const {
  return Point(coordinates_ - translation);
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
constexpr bool Point<Vector>::operator==(Point<Vector> const& right) const {
  return coordinates_ == right.coordinates_;
}

template<typename Vector>
constexpr bool Point<Vector>::operator!=(Point<Vector> const& right) const {
  return coordinates_ != right.coordinates_;
}

template<typename Vector>
void Point<Vector>::WriteToMessage(
    not_null<serialization::Point*> const message) const {
  PointSerializer<Vector>::WriteToMessage(coordinates_, message);
}

template<typename Vector>
template<typename, typename>
Point<Vector> Point<Vector>::ReadFromMessage(
    serialization::Point const& message) {
  Point result;
  result.coordinates_ = PointSerializer<Vector>::ReadFromMessage(message);
  return result;
}

template<typename Vector>
constexpr Point<Vector>::Point(Vector const& coordinates)
    : coordinates_(coordinates) {}

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point) {
  return point + translation;
}

template<typename Vector>
constexpr typename std::enable_if_t<is_quantity_v<Vector>, bool> operator<(
    Point<Vector> const& left,
    Point<Vector> const& right) {
  return left.coordinates_ < right.coordinates_;
}

template<typename Vector>
constexpr typename std::enable_if_t<is_quantity_v<Vector>, bool>
operator<=(Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ <= right.coordinates_;
}

template<typename Vector>
constexpr typename std::enable_if_t<is_quantity_v<Vector>, bool>
operator>=(Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ >= right.coordinates_;
}

template<typename Vector>
constexpr typename std::enable_if_t<is_quantity_v<Vector>, bool> operator>(
    Point<Vector> const& left,
    Point<Vector> const& right) {
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

}  // namespace internal_point

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
  Point<Vector> const origin;
  return origin + weighted_sum_ / weight_;
}

template<typename Vector, typename Weight>
Weight const& BarycentreCalculator<Point<Vector>, Weight>::weight() const {
  return weight_;
}

}  // namespace geometry
}  // namespace principia
