#pragma once

#include <vector>

namespace principia {
namespace geometry {

// Point<Vector> is an affine space on the vector space Vector. Vector should
// be equipped with operators +, -, +=, -=, ==, !=, as well as Vector * Weight
// and Vector / Weight for any Weight used in Barycentre.
template<typename Vector>
class Point {
 public:
  Point() = default;
  explicit Point(Vector const& coordinates);
  ~Point() = default;

  Vector operator-(Point const& from) const;

  Point operator+(Vector const& translation) const;
  Point operator-(Vector const& translation) const;

  Point& operator+=(Vector const& translation);
  Point& operator-=(Vector const& translation);

  bool operator==(Point const& right) const;
  bool operator!=(Point const& right) const;

 private:
  Vector coordinates_;

  template<typename FromFrame, typename ToFrame, typename Scalar,
           template<typename, typename> class LinearMap>
  friend class AffineMap;

  template<typename V>
  friend Point<V> operator+(V const& translation,
                            Point<V> const& point);

  template<typename Vector, typename Weight>
  friend Point<Vector> Barycentre(std::vector<Point<Vector>> const& points,
                                  std::vector<Weight> const& weights);
};

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point);

template<typename Vector, typename Weight>
Point<Vector> Barycentre(std::vector<Point<Vector>> const& points,
                         std::vector<Weight> const& weights);

}  // namespace geometry
}  // namespace principia

#include "geometry/point_body.hpp"
