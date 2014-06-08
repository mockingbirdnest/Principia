#pragma once

#include "physics/n_body_system.hpp"

namespace principia {
namespace benchmarks {

// Simulates the 18 largest solar system bodies (Pluto and all larger bodies)
// over Sir Isaac Newton's life (1643-01-04T00:00:00Z/1727-03-21T00:00:00Z, in
// Julian date JD2321156.5 to JD2351912.5)
// The integration is performed in  the ICRF/J2000.0 reference frame, the origin
// is the barycentre of the solar system, the reference plane is the ecliptic
// mean and mean equinox of J2000.0.
// All data is from the Jet Propulsion Laboratory's HORIZONS system.
void SolarSystem1643To1727();

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/n_body_system_body.hpp"
