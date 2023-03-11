#include <algorithm>
#include <ctime>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "astronomy/frames.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "numerics/apodization.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/frequency_analysis.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"

#if PRINCIPIA_CONTINUOUS_TRAJECTORY_SUPPORTS_PIECEWISE_POISSON_SERIES

namespace principia {
namespace physics {

using namespace principia::astronomy::_frames;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_named_quantities;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::numerics::_fast_fourier_transform;
using namespace principia::numerics::_poisson_series;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::physics::_analytical_series;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

namespace apodization = numerics::apodization;
namespace frequency_analysis = numerics::frequency_analysis;

static constexpr int aperiodic_approximation_degree = 5;
static constexpr int periodic_approximation_degree = 3;
static constexpr int log2_number_of_samples = 14;
static constexpr int number_of_frequencies = 10;
static constexpr Length acceptable_residual = 1 * Metre;
static constexpr Time min_projection_duration = 4 * Day;
static constexpr Time max_projection_duration = 0.25 * JulianYear;

class AnalyticalSeriesTest : public ::testing::Test {
 protected:
  AnalyticalSeriesTest()
      : logger_(TEMP_DIR / "analytical_series.wl",
               /*make_unique=*/false) {
    google::LogToStderr();
  }

  std::string ApplyParameters(std::string const& name) {
    return mathematica::Evaluate(
        mathematica::Apply(name,
                           std::tuple{aperiodic_approximation_degree,
                                      periodic_approximation_degree,
                                      celestial_}));
  }

  template<int degree>
  PoissonSeries<Displacement<ICRS>,
                aperiodic_approximation_degree, periodic_approximation_degree,
                EstrinEvaluator>
  ComputeCompactRepresentation(ContinuousTrajectory<ICRS> const& trajectory) {
    Instant const t_min = trajectory.t_min();
    Instant const t_max = trajectory.t_max();
    auto const piecewise_poisson_series =
        trajectory.ToPiecewisePoissonSeries<degree, 0>(t_min, t_max);

    logger_.Append("tMin",
                   std::tuple(projection_duration_, t_min),
                   mathematica::ExpressIn(Second));
    logger_.Append("tMax",
                   std::tuple(projection_duration_, t_max),
                   mathematica::ExpressIn(Second));

    int step = 0;
    std::clock_t const start_cpu = std::clock();

    auto angular_frequency_calculator =
        [this, start_cpu, &step, t_min, t_max](
            auto const& residual) -> std::optional<AngularFrequency> {
      Time const Δt = (t_max - t_min) / (1 << log2_number_of_samples);
      LOG(INFO) << "step=" << step;
      if (step == 0) {
        ++step;
        return AngularFrequency();
      } else if (step <= number_of_frequencies) {
        ++step;
        Length max_residual;
        std::vector<Displacement<ICRS>> residuals;
        for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
          residuals.push_back(residual(t_min + i * Δt));
          max_residual = std::max(max_residual, residuals.back().Norm());
        }
        LOG(INFO) << "max_residual=" << max_residual;

        std::clock_t const end_cpu = std::clock();
        logger_.Append(
            ApplyParameters("maxResidual"),
            std::tuple(projection_duration_, step - 1,
                       max_residual,
                       1000.0 * (end_cpu - start_cpu) / CLOCKS_PER_SEC),
            mathematica::ExpressIn(Metre, Second));

        if (max_residual < acceptable_residual) {
          return std::nullopt;
        }
        auto fft =
            std::make_unique<FastFourierTransform<Displacement<ICRS>,
                                                  Instant,
                                                  1 << log2_number_of_samples>>(
                residuals, Δt);
        auto const mode = fft->Mode(2 * π * Radian / (t_max - t_min),
                                    Infinity<AngularFrequency>);
        Interval<Time> const period{2 * π * Radian / mode.max,
                                    2 * π * Radian / mode.min};
        LOG(INFO) << "period=" << period;
        auto const precise_mode = frequency_analysis::PreciseMode(
            mode, residual, apodization::Hann<EstrinEvaluator>(t_min, t_max));
        auto const precise_period = 2 * π * Radian / precise_mode;
        LOG(INFO) << "precise_period=" << precise_period;
        return precise_mode;
      } else {
        Length max_residual;
        for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
          max_residual =
              std::max(max_residual, residual(t_min + i * Δt).Norm());
        }
        LOG(INFO) << "max_residual=" << max_residual
                  << (max_residual > acceptable_residual ? " ***" : "");

        std::clock_t const end_cpu = std::clock();
        logger_.Append(
            ApplyParameters("maxResidual"),
            std::tuple(projection_duration_, step,
                       max_residual,
                       1000.0 * (end_cpu - start_cpu) / CLOCKS_PER_SEC),
            mathematica::ExpressIn(Metre, Second));

        return std::nullopt;
      }
    };

    return frequency_analysis::IncrementalProjection<
        aperiodic_approximation_degree, periodic_approximation_degree>(
        piecewise_poisson_series,
        angular_frequency_calculator,
        apodization::Dirichlet<EstrinEvaluator>(t_min, t_max),
        t_min, t_max);
  }

  mathematica::Logger logger_;
  std::string celestial_;
  Time projection_duration_;
};

#define PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(                 \
    degree, approximation, trajectory)                                 \
  case degree: {                                                       \
    approximation =                                                    \
        std::make_unique<PoissonSeries<Displacement<ICRS>,             \
                                       aperiodic_approximation_degree, \
                                       periodic_approximation_degree,  \
                                       EstrinEvaluator>>(              \
            ComputeCompactRepresentation<(degree)>(trajectory));       \
    break;                                                             \
  }

#if !_DEBUG
TEST_F(AnalyticalSeriesTest, CompactRepresentation) {
  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  for (Time projection_duration = min_projection_duration;
       projection_duration <= max_projection_duration;
       projection_duration *= 2) {
    // NOTE(phl): Keep these parameters aligned with
    // sol_numerics_blueprint.proto.txt.
    auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute));
    ephemeris->Prolong(solar_system_at_j2000.epoch() + projection_duration);

    for (auto const& celestial : {"Io", "Moon", "Phobos"}) {
      LOG(INFO) << "---------- " << celestial << " " << projection_duration;
      auto const& celestial_trajectory =
          solar_system_at_j2000.trajectory(*ephemeris, celestial);
      int const celestial_piecewise_poisson_series_degree =
          celestial_trajectory.PiecewisePoissonSeriesDegree(
              celestial_trajectory.t_min(), celestial_trajectory.t_max());
      std::unique_ptr<PoissonSeries<Displacement<ICRS>,
                                    aperiodic_approximation_degree,
                                    periodic_approximation_degree,
                                    EstrinEvaluator>>
          celestial_approximation;

      celestial_ = celestial;
      projection_duration_ = projection_duration;

      switch (celestial_piecewise_poisson_series_degree) {
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            3, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            4, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            5, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            6, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            7, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            8, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            9, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            10, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            11, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            12, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            13, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            14, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            15, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            16, celestial_approximation, celestial_trajectory);
        PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
            17, celestial_approximation, celestial_trajectory);
        default:
          LOG(FATAL) << "Unexpected degree "
                     << celestial_piecewise_poisson_series_degree;
      };

      logger_.Set(ApplyParameters("approximation"),
                  *celestial_approximation,
                  mathematica::ExpressIn(Metre, Second, Radian));
    }
  }
}
#endif

#undef PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE

}  // namespace physics
}  // namespace principia

#endif
