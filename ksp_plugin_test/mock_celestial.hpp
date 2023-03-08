#pragma once

#include "ksp_plugin/celestial.hpp"

#include "gmock/gmock.h"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace _celestial {
namespace internal {

using namespace principia::testing_utilities::_make_not_null;

class MockCelestial : public Celestial {
 public:
  MockCelestial()
      : Celestial(make_not_null<RotatingBody<Barycentric> const*>()) {}

  MOCK_METHOD(DegreesOfFreedom<Barycentric>,
              current_degrees_of_freedom,
              (Instant const& current_time),
              (const, override));
  MOCK_METHOD(Position<Barycentric>,
              current_position,
              (Instant const& current_time),
              (const, override));
  MOCK_METHOD(Velocity<Barycentric>,
              current_velocity,
              (Instant const& current_time),
              (const, override));
};

}  // namespace internal

using internal::MockCelestial;

}  // namespace _celestial
}  // namespace ksp_plugin
}  // namespace principia
