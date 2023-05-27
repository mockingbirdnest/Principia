#pragma once

#include <vector>

#include "physics/equipotential.hpp"

namespace principia {
namespace physics {
namespace _lagrange_equipotentials {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_equipotential;
using namespace principia::physics::_massive_body;

template<typename Inertial, typename RotatingPulsating>
class LagrangeEquipotentials {
 public:
  LagrangeEquipotentials(not_null<Ephemeris<Inertial> const*> ephemeris);

  struct Parameters {
    std::vector<not_null<MassiveBody const*>> primaries;
    std::vector<not_null<MassiveBody const*>> secondaries;
    Instant time;
    // The number of energy levels at which to draw equipotentials.
    int levels = 8;
    // The level of the equipotential going through the L₁ Lagrange point.
    int l1_level = 7;
    // Whether to draw an additional equipotential going through the L₂
    // point.
    bool show_l245_level = true;
  };

  absl::StatusOr<typename Equipotential<Inertial, RotatingPulsating>::Lines>
  ComputeLines(Parameters const& parameters);

  struct Equipotentials {
    Equipotential<Inertial, RotatingPulsating>::Lines lines;
    Parameters parameters;
  };

 private:
  not_null<Ephemeris<Inertial> const*> const ephemeris_;

};

}  // namespace internal

using internal::LagrangeEquipotentials;

}  // namespace _lagrange_equipotentials
}  // namespace physics
}  // namespace principia

#include "physics/lagrange_equipotentials_body.hpp"
