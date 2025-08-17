// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_filter=Geopotential  // NOLINT(whitespace/line_length)

#include <random>
#include <vector>

#include "astronomy/fortran_astrodynamics_toolkit.hpp"  // 🧙 For fortran_astrodynamics_toolkit_.  // NOLINT
#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/legendre_normalization_factor.mathematica.h"
#include "physics/geopotential.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_legendre_normalization_factor;
using namespace principia::physics::_geopotential;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_parser;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

template<typename Frame>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
GeneralSphericalHarmonicsAccelerationCpp(
    Geopotential<Frame> const& geopotential,
    Instant const& t,
    Displacement<Frame> const& r) {
  auto const r² = r.Norm²();
  auto const r_norm = Sqrt(r²);
  auto const one_over_r³ = r_norm / (r² * r²);
  return geopotential.GeneralSphericalHarmonicsAcceleration(
      t, r, r_norm, r², one_over_r³);
}

// For fairness, the Fortran implementation is wrapped to have the same API as
// the C++ one.
template<typename Frame, int degree, int order>
Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
GeneralSphericalHarmonicsAccelerationF90(
    not_null<OblateBody<Frame> const*> const body,
    double const mu,
    double const rbar,
    FixedMatrix<double, degree + 1, order + 1> const& cnm,
    FixedMatrix<double, degree + 1, order + 1> const& snm,
    Instant const& t,
    Displacement<Frame> const& r) {
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;
  auto const from_surface_frame =
      body->template FromSurfaceFrame<SurfaceFrame>(t);
  auto const to_surface_frame = from_surface_frame.Inverse();

  Displacement<SurfaceFrame> const r_surface = to_surface_frame(r);
  auto const acceleration_surface =
      Vector<Quotient<Acceleration, GravitationalParameter>, SurfaceFrame>(
          si::Unit<Quotient<Acceleration, GravitationalParameter>> *
          astronomy::_fortran_astrodynamics_toolkit::
              ComputeGravityAccelerationLear<degree, order>(
                  r_surface.coordinates() / Metre, mu, rbar, cnm, snm));
  return from_surface_frame(acceleration_surface);
}

OblateBody<ICRS> MakeEarthBody(SolarSystem<ICRS>& solar_system,
                               int const max_degree) {
  solar_system.LimitOblatenessToDegree("Earth", max_degree);
  auto earth_message = solar_system.gravity_model_message("Earth");

  Angle const earth_right_ascension_of_pole = 0 * Degree;
  Angle const earth_declination_of_pole = 90 * Degree;
  auto const earth_μ = solar_system.gravitational_parameter("Earth");
  auto const earth_reference_radius =
      ParseQuantity<Length>(earth_message.reference_radius());
  MassiveBody::Parameters const massive_body_parameters(earth_μ);
  RotatingBody<ICRS>::Parameters rotating_body_parameters(
      /*mean_radius=*/solar_system.mean_radius("Earth"),
      /*reference_angle=*/0 * Radian,
      /*reference_instant=*/Instant(),
      /*angular_frequency=*/1 * Radian / Second,
      earth_right_ascension_of_pole,
      earth_declination_of_pole);
  return OblateBody<ICRS>(
      massive_body_parameters,
      rotating_body_parameters,
      OblateBody<ICRS>::Parameters::ReadFromMessage(
          earth_message.geopotential(), earth_reference_radius));
}

void BM_ComputeGeopotentialCpp(benchmark::State& state) {
  int const max_degree = state.range(0);

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");

  auto const earth = MakeEarthBody(solar_system_2000, max_degree);
  Geopotential<ICRS> const geopotential(&earth, /*tolerance=*/0);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution(-1e7, 1e7);
  std::vector<Displacement<ICRS>> displacements;
  for (int i = 0; i < 1e3; ++i) {
    displacements.push_back(earth.FromSurfaceFrame<ITRS>(Instant())(
        Displacement<ITRS>({distribution(random) * Metre,
                            distribution(random) * Metre,
                            distribution(random) * Metre})));
  }

  for (auto _ : state) {
    Vector<Exponentiation<Length, -2>, ICRS> acceleration;
    for (auto const& displacement : displacements) {
      acceleration = GeneralSphericalHarmonicsAccelerationCpp(
                         geopotential, Instant(), displacement);
    }
    benchmark::DoNotOptimize(acceleration);
  }
}

void BM_ComputeGeopotentialDistance(benchmark::State& state) {
  // Check the performance around this distance.  May be used to tell apart the
  // various contributions.
  double const distance_in_kilometres = state.range(0);

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");

  auto const earth = MakeEarthBody(solar_system_2000, /*max_degree=*/10);
  Geopotential<ICRS> const geopotential(&earth, /*tolerance=*/0x1.0p-24);

  // Generate points in a spherical shell.
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution(0 * distance_in_kilometres,
                                                1.1 * distance_in_kilometres);
  std::vector<Displacement<ICRS>> displacements;
  for (int i = 0; i < 1e3;) {
    Displacement<ITRS> const displacement({distribution(random) * Kilo(Metre),
                                           distribution(random) * Kilo(Metre),
                                           distribution(random) * Kilo(Metre)});
    if (displacement.Norm() > 0.9 * distance_in_kilometres * Kilo(Metre) &&
        displacement.Norm() < 1.1 * distance_in_kilometres * Kilo(Metre)) {
      displacements.push_back(
          earth.FromSurfaceFrame<ITRS>(Instant())(displacement));
      ++i;
    }
  }

  for (auto _ : state) {
    Vector<Exponentiation<Length, -2>, ICRS> acceleration;
    for (auto const& displacement : displacements) {
      acceleration = GeneralSphericalHarmonicsAccelerationCpp(
                         geopotential, Instant(), displacement);
    }
    benchmark::DoNotOptimize(acceleration);
  }
}

#define PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(d)                         \
  case (d): {                                                              \
    FixedMatrix<double, (d) + 1, (d) + 1> cnm;                             \
    FixedMatrix<double, (d) + 1, (d) + 1> snm;                             \
    for (int n = 0; n <= (d); ++n) {                                       \
      for (int m = 0; m <= n; ++m) {                                       \
        cnm(n, m) = earth.cos()(n, m) * LegendreNormalizationFactor(n, m); \
        snm(n, m) = earth.sin()(n, m) * LegendreNormalizationFactor(n, m); \
      }                                                                    \
    }                                                                      \
    for (auto _ : state) {                                                 \
      Vector<Exponentiation<Length, -2>, ICRS> acceleration;               \
      for (auto const& displacement : displacements) {                     \
        acceleration =                                                     \
            GeneralSphericalHarmonicsAccelerationF90<ICRS, (d), (d)>(      \
                &earth, mu, rbar, cnm, snm, Instant(), displacement);      \
      }                                                                    \
      benchmark::DoNotOptimize(acceleration);                              \
    }                                                                      \
    break;                                                                 \
  }

void BM_ComputeGeopotentialF90(benchmark::State& state) {
  int const max_degree = state.range(0);

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const earth = MakeEarthBody(solar_system_2000, max_degree);

  double mu =
      earth.gravitational_parameter() / si::Unit<GravitationalParameter>;
  double rbar = earth.reference_radius() / Metre;

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution(-1e7, 1e7);
  std::vector<Displacement<ICRS>> displacements;
  for (int i = 0; i < 1e3; ++i) {
    displacements.push_back(earth.FromSurfaceFrame<ITRS>(Instant())(
        Displacement<ITRS>({distribution(random) * Metre,
                            distribution(random) * Metre,
                            distribution(random) * Metre})));
  }

  switch (max_degree) {
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(2);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(3);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(4);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(5);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(6);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(7);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(8);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(9);
    PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90(10);
  }
}

#undef PRINCIPIA_CASE_COMPUTE_GEOPOTENTIAL_F90

BENCHMARK(BM_ComputeGeopotentialCpp)
    ->Arg(2)
    ->Arg(3)
    ->Arg(5)
    ->Arg(10)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_ComputeGeopotentialF90)
    ->Arg(2)
    ->Arg(3)
    ->Arg(5)
    ->Arg(10)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_ComputeGeopotentialDistance)
    ->Arg(150'000)    // C₂₂, S₂₂, J₂.
    ->Arg(500'000)    // J₂.
    ->Arg(5'000'000)  // Central
    ->Unit(benchmark::kMicrosecond);

}  // namespace physics
}  // namespace principia
