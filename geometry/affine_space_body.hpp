#pragma once

namespace principia {
namespace geometry {

template<typename Vector>
Point<Vector>::Point(Vector const& coordinates) : coordinates_(coordinates) {};

template<typename Vector>
Vector operator-(Point<Vector> const& to, Point<Vector> const& from) {
  return to.coordinates_ - from.coordinates_;
}

template<typename Vector>
Point<Vector> operator+(Point<Vector> const& left, Vector const& right) {
  return Point<Vector>(left.coordinates_ + right);
}

template<typename Vector>
Point<Vector> operator+(Vector const& left, Point<Vector> const& right) {
  return Point<Vector>(left + right.coordinates_);
}

template<typename Vector>
Point<Vector> operator-(Point<Vector> const& left, Vector const& right) {
  return Point<Vector>(left.coordinates_ - right);
}

template<typename Vector>
bool operator==(Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ == right.coordinates_;
}

template<typename Vector>
bool operator!=(Point<Vector> const& left, Point<Vector> const& right) {
  return left.coordinates_ != right.coordinates_;
}

template<typename Vector, typename Weight>
Point<Vector> Barycentre(Point<Vector> const& left,
                         Weight const& left_weight,
                         Point<Vector> const& right,
                         Weight const& right_weight) {
  return Point<Vector>(
      (left.coordinates_ * left_weight + right.coordinates_ * right_weight) /
          (left_weight + right_weight));
}

}  // namespace geometry
}  // namespace principia
