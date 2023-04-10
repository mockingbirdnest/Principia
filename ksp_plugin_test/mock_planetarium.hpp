#pragma once

#include "ksp_plugin/planetarium.hpp"

#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace _planetarium {
namespace internal {

using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_make_not_null;

class MockPlanetarium : public Planetarium {
 public:
  MockPlanetarium()
      : Planetarium(
            Planetarium::Parameters(1.0, 1.0 * Radian, 1.0 * Radian),
            Perspective<Navigation, Camera>(
                RigidTransformation<Navigation, Camera>(
                    Navigation::origin,
                    Camera::origin,
                    Signature<Navigation, Camera>::CentralInversion()
                    .Forget<OrthogonalMap>()).Forget<Similarity>(),
                1 * Metre),
            make_not_null<Ephemeris<Barycentric> const*>(),
            make_not_null<NavigationFrame const*>(),
            [](Position<Navigation> const& plotted_point) {
              constexpr auto inverse_scale_factor = 1 / (6000 * Metre);
              return ScaledSpacePoint::FromCoordinates(
                  ((plotted_point - Navigation::origin) *
                   inverse_scale_factor).coordinates());
            }) {}
};

}  // namespace internal

using internal::MockPlanetarium;

}  // namespace _planetarium
}  // namespace ksp_plugin
}  // namespace principia
