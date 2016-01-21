#include "jacobi_coordinates.hpp"

#include "astronomy/frames.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;

namespace physics {

template class JacobiCoordinates<ICRFJ2000Equator>;
template class HierarchicalSystem<ICRFJ2000Equator>;

}  // namespace physics
}  // namespace principia
