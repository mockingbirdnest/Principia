#pragma once

#include "ksp_plugin/celestial.hpp"

#include "gmock/gmock.h"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_celestial {

using testing_utilities::make_not_null;

class MockCelestial : public Celestial {
 public:
  MockCelestial()
      : Celestial(make_not_null<RotatingBody<Barycentric> const*>()) {}

  MOCK_CONST_METHOD1(
      current_degrees_of_freedom,
      DegreesOfFreedom<Barycentric>(Instant const& current_time));
  MOCK_CONST_METHOD1(current_position,
                     Position<Barycentric>(Instant const& current_time));
  MOCK_CONST_METHOD1(current_velocity,
                     Velocity<Barycentric>(Instant const& current_time));
};

}  // namespace internal_celestial

using internal_celestial::MockCelestial;

}  // namespace ksp_plugin
}  // namespace principia
