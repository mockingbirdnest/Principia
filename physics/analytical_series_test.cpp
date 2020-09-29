
#include "astronomy/frames.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/apodization.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/frequency_analysis.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"

namespace principia {
namespace physics {

using astronomy::ICRS;
using geometry::Displacement;
using geometry::Instant;
using geometry::Interval;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using numerics::EstrinEvaluator;
using numerics::FastFourierTransform;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace apodization = numerics::apodization;
namespace frequency_analysis = numerics::frequency_analysis;

class AnalyticalSeriesTest : public ::testing::Test {
 protected:
  AnalyticalSeriesTest() {
    google::LogToStderr();
  }
};

TEST_F(AnalyticalSeriesTest, SolarSystemSeries) {
  mathematica::Logger logger(TEMP_DIR / "analytical_series.wl",
                             /*make_unique=*/false);
  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  // NOTE(phl): Keep these parameters aligned with
  // sol_numerics_blueprint.proto.txt.
  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
          /*step=*/10 * Minute));
  ephemeris->Prolong(solar_system_at_j2000.epoch() + 0.25 * JulianYear);

  auto const& io_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Io");
  Instant const t_min = io_trajectory.t_min();
  Instant const t_max = io_trajectory.t_max();
  int const io_piecewise_poisson_series_degree =
      io_trajectory.PiecewisePoissonSeriesDegree(t_min, t_max);
  LOG(ERROR)<<io_piecewise_poisson_series_degree;

  // TODO(phl): Use a switch statement with macros.
  auto const io_piecewise_poisson_series =
      io_trajectory.ToPiecewisePoissonSeries<7>(t_min, t_max);
  logger.Set("tMin", t_min, mathematica::ExpressIn(Second));
  logger.Set("tMax", t_max, mathematica::ExpressIn(Second));

  std::vector<Displacement<ICRS>> displacements;
  std::vector<std::tuple<Instant, Displacement<ICRS>>> trajectory;
  for (int i = 0; i <= 1000; ++i) {
    auto const t = t_min + i * (t_max - t_min) / 1000;
    auto const current_displacements = io_piecewise_poisson_series(t);
    displacements.push_back(current_displacements);
    auto const current_trajectory =
        io_trajectory.EvaluatePosition(t) - ICRS::origin;
    trajectory.push_back({t, current_trajectory});
  }

  static constexpr int number_of_frequencies = 10;
  static constexpr int log2_number_of_samples = 10;

  int step = 0;
  auto angular_frequency_calculator =
      [&logger, &step, t_min, t_max](
          auto const& residual) -> std::optional<AngularFrequency> {
    LOG(ERROR) << "step=" << step;
    if (step == 0) {
      ++step;
      return AngularFrequency();
    } else if (step <= number_of_frequencies) {
      ++step;
      Length max_residual;
      std::vector<Displacement<ICRS>> residuals;
      Time const Δt = (t_max - t_min) / (1 << log2_number_of_samples);
      for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
        residuals.push_back(residual(t_min + i * Δt));
        max_residual = std::max(max_residual, residuals.back().Norm());
      }
      LOG(ERROR) << "max_residual=" << max_residual;
      auto fft =
          std::make_unique<FastFourierTransform<Displacement<ICRS>,
                                                1 << log2_number_of_samples>>(
              residuals, Δt);
      auto const mode = fft->Mode();
      Interval<Time> const period{2 * π * Radian / mode.max,
                                  2 * π * Radian / mode.min};
      LOG(ERROR) << "period=" << period;
      auto const precise_mode =
          PreciseMode(mode,
                      residual,
                      apodization::Hann<EstrinEvaluator>(t_min, t_max));
      auto const precise_period = 2 * π * Radian / precise_mode;
      LOG(ERROR) << "precise_period=" << precise_period;
      logger.Append(
          "precisePeriods", precise_period, mathematica::ExpressIn(Second));
      return precise_mode;
    } else {
      Length max_residual;
      std::vector<Displacement<ICRS>> residuals;
      Time const Δt = (t_max - t_min) / (1 << log2_number_of_samples);
      for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
        residuals.push_back(residual(t_min + i * Δt));
        max_residual = std::max(max_residual, residuals.back().Norm());
      }
      LOG(ERROR) << "max_residual=" << max_residual;
      return std::nullopt;
    }
  };

  auto const solution = frequency_analysis::IncrementalProjection<5>(
      io_piecewise_poisson_series,
      angular_frequency_calculator,
      apodization::Dirichlet<EstrinEvaluator>(t_min, t_max),
      t_min, t_max);
}

}  // namespace physics
}  // namespace principia
