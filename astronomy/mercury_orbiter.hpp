#pragma once

#include "astronomy/time_scales.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace astronomy {
namespace internal_mercury_orbiter {

using astronomy::operator""_TT;
using physics::DegreesOfFreedom;
using quantities::si::Metre;
using quantities::si::Second;
using namespace principia::geometry::_named_quantities;

// State of the spacecraft Mercury Orbiter 1 at the start of September, from a
// save by Butcher given in issue #1119.  This is used both from a plugin
// compatibility test which reads that save, and from an astronomical test which
// looks for the Лидов–古在 mechanism in the orbit.

constexpr Instant MercuryOrbiterInitialTime =
    "1966-09-01T00:16:55"_TT + 0.2571494579315186 * Second;

template<typename Barycentric>
physics::DegreesOfFreedom<Barycentric> const
    MercuryOrbiterInitialDegreesOfFreedom = {
        Barycentric::origin +
            Displacement<Barycentric>({-2.40627773705000000e+10 * Metre,
                                       +3.52445087251250000e+10 * Metre,
                                       +2.13640458684375000e+10 * Metre}),
        Velocity<Barycentric>({-5.19594188203811646e+04 * (Metre / Second),
                               -2.23741500134468079e+04 * (Metre / Second),
                               -7.15344990825653076e+03 * (Metre / Second)})};

}  // namespace internal_mercury_orbiter

using internal_mercury_orbiter::MercuryOrbiterInitialDegreesOfFreedom;
using internal_mercury_orbiter::MercuryOrbiterInitialTime;

}  // namespace astronomy
}  // namespace principia
