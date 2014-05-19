#pragma once

namespace principia {
namespace geometry {

template<typename Vector>
Point<Vector>::Point(Vector const& coordinates) : coordinates_(coordinates) {};

template<typename Vector>
Vector Point<Vector>::operator-(Point<Vector> const& from) const {
  return coordinates_ - from.coordinates_;
}

template<typename Vector>
Point<Vector> Point<Vector>::operator+(Vector const& translation) const {
  return Point<Vector>(coordinates_ + translation);
}

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point) {
  return Point<Vector>(translation + point.coordinates_);
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
