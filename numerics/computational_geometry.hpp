#include "geometry/point.hpp"

namespace principia {
namespace numerics {
namespace internal_computational_geometry {

using geometry::Point;

template<typename Vector>
class NearestNeighbourCalculator {
 public:
  NearestNeighbourCalculator(int points_per_cell,
    int initial_capacity,
                             std::pair<Point<Vector>> const& corners);

  void Insert(Point<Vector> const& point);

  Point<Vector> const& Search(Point<Vector> const& point) const;
};

}  // namespace internal_computational_geometry
}  // namespace numerics
}  // namespace principia

#include "numerics/computational_geometry_body.hpp"
