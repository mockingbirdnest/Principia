#include "numerics/computational_geometry.hpp"

namespace principia {
namespace numerics {
namespace internal_computational_geometry {

template<typename Vector>
NearestNeighbourCalculator<Vector>::NearestNeighbourCalculator() {}

template<typename Vector>
void NearestNeighbourCalculator<Vector>::Insert(Point<Vector> const& point) {}

template<typename Vector>
Point<Vector> const& NearestNeighbourCalculator<Vector>::Search(
    Point<Vector> const& point) const {
  // TODO: insert return statement here
}

}  // namespace internal_computational_geometry
}  // namespace numerics
}  // namespace principia
