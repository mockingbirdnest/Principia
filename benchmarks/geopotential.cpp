
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Geopotential  // NOLINT(whitespace/line_length)

#include "physics/geopotential_body.hpp"

#include <random>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using astronomy::ICRS;
using astronomy::ITRS;
using geometry::Displacement;
using geometry::Instant;
using geometry::Vector;
using physics::SolarSystem;
using quantities::Acceleration;
using quantities::Angle;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::ParseQuantity;
using quantities::Pow;
using quantities::Quotient;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
GeneralSphericalHarmonicsAcceleration(Geopotential<Frame> const& geopotential,
                                      Instant const& t,
                                      Displacement<Frame> const& r) {
  auto const r² = r.Norm²();
  auto const one_over_r³ = 1.0 / (r² * r.Norm());
  return geopotential.GeneralSphericalHarmonicsAcceleration(
              t, r, r², one_over_r³);
}

void BM_ComputeGeopotential(benchmark::State& state) {
  int const max_degree = state.range_x();

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto earth_message = solar_system_2000.gravity_model_message("Earth");
  earth_message.mutable_geopotential()->set_max_degree(max_degree);

  Angle const earth_right_ascension_of_pole = 0 * Degree;
  Angle const earth_declination_of_pole = 90 * Degree;
  auto const earth_μ = solar_system_2000.gravitational_parameter("Earth");
  auto const earth_reference_radius =
      ParseQuantity<Length>(earth_message.reference_radius());
  MassiveBody::Parameters const massive_body_parameters(earth_μ);
  RotatingBody<ICRS>::Parameters rotating_body_parameters(
      /*mean_radius=*/solar_system_2000.mean_radius("Earth"),
      /*reference_angle=*/0 * Radian,
      /*reference_instant=*/Instant(),
      /*angular_frequency=*/1 * Radian / Second,
      earth_right_ascension_of_pole,
      earth_declination_of_pole);
  OblateBody<ICRS> const earth = OblateBody<ICRS>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
  Geopotential<ICRS> const geopotential(&earth);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> const distribution(-1e7, 1e7);
  std::vector<Displacement<ICRS>> displacements;
  for (int i = 0; i < 1e3; ++i) {
    displacements.push_back(earth.FromSurfaceFrame<ITRS>(Instant())(
        Displacement<ITRS>({distribution(random) * Metre,
                            distribution(random) * Metre,
                            distribution(random) * Metre})));
  }

  while (state.KeepRunning()) {
    Vector<Acceleration, ICRS> acceleration;
    for (auto const& displacement : displacements) {
      acceleration = earth_μ * GeneralSphericalHarmonicsAcceleration(
                                   geopotential, Instant(), displacement);
    }
    benchmark::DoNotOptimize(acceleration);
  }
}

BENCHMARK(BM_ComputeGeopotential)->Arg(2)->Arg(3)->Arg(5)->Arg(10);

}  // namespace physics
}  // namespace principia
