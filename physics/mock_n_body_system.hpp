#pragma once

#include "physics/ephemeris.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace physics {

template<typename InertialFrame>
class MockNBodySystem : public Ephemeris<InertialFrame> {
 public:
  MockNBodySystem() : Ephemeris() {}

  // TODO(egg): MOCK_ALL_THE_THINGS
};

}  // namespace physics
}  // namespace principia
