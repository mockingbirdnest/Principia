
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Geopotential  // NOLINT(whitespace/line_length)

#include "physics/geopotential_body.hpp"

#include <random>
#include <vector>

#include "astronomy/fortran_astrodynamics_toolkit.hpp"
#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/legendre.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {

using astronomy::ICRS;
using astronomy::ITRS;
using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::R3Element;
using geometry::Vector;
using numerics::FixedMatrix;
using numerics::LegendreNormalizationFactor;
using physics::SolarSystem;
using quantities::Acceleration;
using quantities::Angle;
using quantities::Exponentiation;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::ParseQuantity;
using quantities::Pow;
using quantities::Quotient;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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
  struct SurfaceFrame;
  auto const from_surface_frame =
      body->template FromSurfaceFrame<SurfaceFrame>(t);
  auto const to_surface_frame = from_surface_frame.Inverse();

  Displacement<SurfaceFrame> const r_surface = to_surface_frame(r);
  auto const acceleration_surface =
      Vector<Quotient<Acceleration, GravitationalParameter>, SurfaceFrame>(
          SIUnit<Quotient<Acceleration, GravitationalParameter>>() *
          astronomy::fortran_astrodynamics_toolkit::
              ComputeGravityAccelerationLear<degree, order>(
                  r_surface.coordinates() / Metre, mu, rbar, cnm, snm));
  return from_surface_frame(acceleration_surface);
}

OblateBody<ICRS> MakeEarthBody(SolarSystem<ICRS> const& solar_system,
                               int const max_degree) {
  auto earth_message = solar_system.gravity_model_message("Earth");
  earth_message.mutable_geopotential()->set_max_degree(max_degree);
  earth_message.mutable_geopotential()->clear_zonal();

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
  int const max_degree = state.range_x();

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");

  auto const earth = MakeEarthBody(solar_system_2000, max_degree);
  Geopotential<ICRS> const geopotential(&earth, /*tolerance=0*/);

  std::mt19937_64 random(42);
  std::uniform_real_distribution<> distribution(-1e7, 1e7);
  std::vector<Displacement<ICRS>> displacements;
  for (int i = 0; i < 1e3; ++i) {
    displacements.push_back(earth.FromSurfaceFrame<ITRS>(Instant())(
        Displacement<ITRS>({distribution(random) * Metre,
                            distribution(random) * Metre,
                            distribution(random) * Metre})));
  }

  while (state.KeepRunning()) {
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
    numerics::FixedMatrix<double, (d) + 1, (d) + 1> cnm;                   \
    numerics::FixedMatrix<double, (d) + 1, (d) + 1> snm;                   \
    for (int n = 0; n <= (d); ++n) {                                       \
      for (int m = 0; m <= n; ++m) {                                       \
        cnm[n][m] = earth.cos()[n][m] * LegendreNormalizationFactor(n, m); \
        snm[n][m] = earth.sin()[n][m] * LegendreNormalizationFactor(n, m); \
      }                                                                    \
    }                                                                      \
    while (state.KeepRunning()) {                                          \
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
  int const max_degree = state.range_x();

  SolarSystem<ICRS> solar_system_2000(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const earth = MakeEarthBody(solar_system_2000, max_degree);

  double mu =
      earth.gravitational_parameter() / SIUnit<GravitationalParameter>();
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

BENCHMARK(BM_ComputeGeopotentialCpp)->Arg(2)->Arg(3)->Arg(5)->Arg(10);
BENCHMARK(BM_ComputeGeopotentialF90)->Arg(2)->Arg(3)->Arg(5)->Arg(10);

}  // namespace physics
}  // namespace principia
