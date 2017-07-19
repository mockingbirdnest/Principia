#pragma once

#include "ksp_plugin/planetarium.hpp"

#include "geometry/affine_map.hpp"
#include "gmock/gmock.h"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_planetarium {

using geometry::AffineMap;
using quantities::si::Metre;
using testing_utilities::make_not_null;

class MockPlanetarium : public Planetarium {
 public:
  MockPlanetarium()
      : Planetarium(Planetarium::Parameters(1.0),
                    Perspective<Navigation, Camera, Length, OrthogonalMap>(
                        AffineMap<Navigation, Camera, Length, OrthogonalMap>::
                            Identity(),
                        1 * Metre),
                    make_not_null<Ephemeris<Barycentric> const*>(),
                    make_not_null<NavigationFrame const*>()) {}
};

}  // namespace internal_planetarium

using internal_planetarium::MockPlanetarium;

}  // namespace ksp_plugin
}  // namespace principia
