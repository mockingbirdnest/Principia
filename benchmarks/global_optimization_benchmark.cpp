// .\Release\x64\benchmarks.exe --benchmark_filter=MLSL --benchmark_repetitions=1  // NOLINT(whitespace/line_length)

#include <cstdint>
#include <memory>

#include "absl/log/globals.h"
#include "absl/log/initialize.h"
#include "absl/strings/str_cat.h"
#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/algebra.hpp"
#include "benchmark/benchmark.h"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/global_optimization.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_view.hpp"
#include "physics/ephemeris.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "physics/trajectory.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/optimization_test_functions.hpp"

namespace principia {
namespace numerics {
namespace {

using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_standard_product_3;
using namespace principia::base::_algebra;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_global_optimization;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_view;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_solar_system;
using namespace principia::physics::_trajectory;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_optimization_test_functions;

using World = Frame<struct WorldTag>;

class GlobalOptimizationBenchmark : public benchmark::Fixture {
 protected:
  // Benchmark doesn't have that capability, so we have to do it ourselves.
  // This function takes about 40 s to run.
  static void SetUpFixture() {
    solar_system_2010_ = std::make_unique<SolarSystem<ICRS>>(
        SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" /
            "sol_initial_state_jd_2455200_500000000.proto.txt").release();
    ephemeris_ =
        solar_system_2010_
            ->MakeEphemeris(
                /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                         /*geopotential_tolerance=*/0x1p-24},
                Ephemeris<ICRS>::FixedStepParameters(
                    SymmetricLinearMultistepIntegrator<
                        QuinlanTremaine1990Order12,
                        Ephemeris<ICRS>::NewtonianMotionEquation>(),
                    /*step=*/10 * Minute))
            .release();
    earth_ = dynamic_cast_not_null<OblateBody<ICRS> const*>(
        solar_system_2010_->massive_body(*ephemeris_, "Earth"));
    earth_trajectory_ = ephemeris_->trajectory(earth_);

    StandardProduct3 const ilrsa_lageos2_sp3(
        SOLUTION_DIR / "astronomy" / "standard_product_3" /
            "ilrsa.orb.lageos2.160319.v35.sp3",
        StandardProduct3::Dialect::ILRSA);
    CHECK(ilrsa_lageos2_sp3.file_has_velocities());
    StandardProduct3::SatelliteIdentifier const lageos2_id{
        StandardProduct3::SatelliteGroup::General, 52};

    auto const& ilrsa_lageos2_trajectory_itrs =
        ilrsa_lageos2_sp3.orbit(lageos2_id).front();
    auto const begin = ilrsa_lageos2_trajectory_itrs->begin();
    CHECK_OK(ephemeris_->Prolong(begin->time));

    BodySurfaceReferenceFrame<ICRS, ITRS> const itrs(ephemeris_, earth_);
    ilrsa_lageos2_trajectory_icrs_ = new DiscreteTrajectory<ICRS>;
    CHECK_OK(ilrsa_lageos2_trajectory_icrs_->Append(
        begin->time,
        itrs.FromThisFrameAtTime(begin->time)(begin->degrees_of_freedom)));
    CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
        ilrsa_lageos2_trajectory_icrs_,
        Ephemeris<ICRS>::NoIntrinsicAcceleration,
        begin->time + 1 * JulianYear,
        Ephemeris<ICRS>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));
  }

  void SetUp(benchmark::State&) override {
    static int const set_up_fixture = []() {
      absl::InitializeLog();
      absl::SetVLogLevel("global_optimization_body", 1);
      absl::SetMinLogLevel(absl::LogSeverityAtLeast::kInfo);
      absl::SetStderrThreshold(absl::LogSeverityAtLeast::kInfo);
      SetUpFixture();
      return 0;
    }();
  }

  // TODO(phl): This code is wrapping an `Instant` into a `Vector` of `Time`,
  // which is ugly as hell, but it will do for now.
  template<typename Frame>
  using Argument = Vector<Time, Frame>;

  template<typename Frame>
  void ComputeApsidesUsingMLSL(Trajectory<Frame> const& reference,
                               DiscreteTrajectoryView<Frame> const& trajectory,
                               int const max_points,
                               std::vector<Argument<Frame>>& apoapsides,
                               std::vector<Argument<Frame>>& periapsides) {
    static Instant const t0;
    Instant const t_min = std::max(reference.t_min(), trajectory.t_min());
    Instant const t_max = std::min(reference.t_max(), trajectory.t_max());
    int f_evaluations = 0;
    int grad_f_evaluations = 0;

    auto const f = [&](Argument<Frame> const& argument) {
      ++f_evaluations;
      Instant const t =
          std::min(std::max(t0 + argument.coordinates()[0], t_min), t_max);
      Position<Frame> const reference_position = reference.EvaluatePosition(t);
      Position<Frame> const trajectory_position =
          trajectory.EvaluatePosition(t);
      Displacement<Frame> const displacement =
          trajectory_position - reference_position;
      return displacement.Norm²();
    };

    auto const grad_f = [&](Argument<Frame> const& argument) {
      ++grad_f_evaluations;
      Instant const t = t0 + argument.coordinates()[0];
      if (t < t_min || t > t_max) {
        return Vector<Variation<Square<Length>>, Frame>();
      }
      DegreesOfFreedom<Frame> const reference_dof =
          reference.EvaluateDegreesOfFreedom(t);
      DegreesOfFreedom<Frame> const trajectory_dof =
          trajectory.EvaluateDegreesOfFreedom(t);
      Displacement<Frame> const displacement =
          trajectory_dof.position() - reference_dof.position();
      return Vector<Variation<Square<Length>>, Frame>(
          {2 * InnerProduct(
                   displacement,
                   trajectory_dof.velocity() - reference_dof.velocity()),
           0 * Pow<2>(Metre) / Second,
           0 * Pow<2>(Metre) / Second});
    };

    using MLSL = MultiLevelSingleLinkage<Square<Length>,
                                         Argument<Frame>,
                                         /*dimensions=*/1>;
    Instant const centre = Barycentre({t_min, t_max}, {1, 1});
    Time const radius = t_max - centre;
    typename MLSL::Box const box{
        Argument<Frame>({centre - t0, 0 * Second, 0 * Second}),
        Argument<Frame>({radius, 0 * Second, 0 * Second})};
    MLSL multi_level_single_linkage(box, f, grad_f);

    periapsides = multi_level_single_linkage.FindGlobalMinima(
        /*points_per_round=*/20,
        /*number_of_rounds=*/std::nullopt,
        /*local_search_tolerance=*/10 * Second);
    apoapsides = multi_level_single_linkage.FindGlobalMaxima(
        /*points_per_round=*/20,
        /*number_of_rounds=*/std::nullopt,
        /*local_search_tolerance=*/10 * Second);
  }

  static SolarSystem<ICRS>* solar_system_2010_;
  static Ephemeris<ICRS>* ephemeris_;
  static OblateBody<ICRS> const* earth_;
  static ContinuousTrajectory<ICRS> const* earth_trajectory_;
  static DiscreteTrajectory<ICRS>* ilrsa_lageos2_trajectory_icrs_;
};

SolarSystem<ICRS>* GlobalOptimizationBenchmark::solar_system_2010_ = nullptr;
Ephemeris<ICRS>* GlobalOptimizationBenchmark::ephemeris_ = nullptr;
OblateBody<ICRS> const* GlobalOptimizationBenchmark::earth_ = nullptr;
ContinuousTrajectory<ICRS> const*
    GlobalOptimizationBenchmark::earth_trajectory_ = nullptr;
DiscreteTrajectory<ICRS>*
    GlobalOptimizationBenchmark::ilrsa_lageos2_trajectory_icrs_ = nullptr;


BENCHMARK_F(GlobalOptimizationBenchmark,
            ComputeApsidesUsingMLSL)(benchmark::State& state) {
  for (auto _ : state) {
    std::vector<Argument<ICRS>> apoapsides;
    std::vector<Argument<ICRS>> periapsides;
    ComputeApsidesUsingMLSL(
        *earth_trajectory_,
        DiscreteTrajectoryView(ilrsa_lageos2_trajectory_icrs_),
        /*max_points=*/std::numeric_limits<int>::max(),
        apoapsides,
        periapsides);
    // TODO(phl): 2 is obviously wrong, somehow the search for periapsides stops
    // too early.  2366 differs from what `ComputeApsides` finds, which is 2364.
    CHECK_EQ(2366, apoapsides.size());
    CHECK_EQ(2, periapsides.size());
  }
}


void BM_MLSLBranin(benchmark::State& state) {
  absl::log_internal::ScopedMinLogLevel s(absl::LogSeverityAtLeast::kError);
  std::int64_t const points_per_round = state.range(0);
  std::int64_t const number_of_rounds = state.range(1);

  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;

  auto branin = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Branin(x₁, x₂);
  };

  auto grad_branin = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁Branin(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0 * Metre, 2.5 * Metre, 7.5 * Metre}),
      .vertices = {
          Displacement<World>({0 * Metre, 7.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 7.5 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, branin, grad_branin);

  int64_t total_minima = 0;
  for (auto _ : state) {
    total_minima +=
        optimizer.FindGlobalMinima(points_per_round,
                                   number_of_rounds, tolerance).size();
  }
  state.SetLabel(
      absl::StrCat("number of minima: ",
                   static_cast<double>(total_minima) / state.iterations()));
}

void BM_MLSLGoldsteinPrice(benchmark::State& state) {
  absl::log_internal::ScopedMinLogLevel s(absl::LogSeverityAtLeast::kError);
  std::int64_t const points_per_round = state.range(0);
  std::int64_t const number_of_rounds = state.range(1);

  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;

  auto goldstein_price = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return GoldsteinPrice(x₁, x₂);
  };

  auto grad_goldstein_price = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁GoldsteinPrice(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>(),
      .vertices = {
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, goldstein_price, grad_goldstein_price);

  int64_t total_minima = 0;
  for (auto _ : state) {
    total_minima +=
        optimizer.FindGlobalMinima(points_per_round,
                                   number_of_rounds, tolerance).size();
  }
  state.SetLabel(
      absl::StrCat("number of minima: ",
                   static_cast<double>(total_minima) / state.iterations()));
}

void BM_MLSLHartmann3(benchmark::State& state) {
  absl::log_internal::ScopedMinLogLevel s(absl::LogSeverityAtLeast::kError);
  std::int64_t const points_per_round = state.range(0);
  std::int64_t const number_of_rounds = state.range(1);

  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/3>;

  auto hartmann3 = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Hartmann3(x₀, x₁, x₂);
  };

  auto grad_hartmann3 = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₀ = coordinates[0] / Metre;
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    auto const [g₀, g₁, g₂] = 𝛁Hartmann3(x₀, x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.5 * Metre}),
      .vertices = {
          Displacement<World>({0.5 * Metre, 0 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 0.5 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, hartmann3, grad_hartmann3);

  int64_t total_minima = 0;
  for (auto _ : state) {
    total_minima +=
        optimizer.FindGlobalMinima(points_per_round,
                                   number_of_rounds, tolerance).size();
  }
  state.SetLabel(
      absl::StrCat("number of minima: ",
                   static_cast<double>(total_minima) / state.iterations()));
}

BENCHMARK(BM_MLSLBranin)->ArgsProduct({{10, 20, 50}, {10, 20, 50}});
BENCHMARK(BM_MLSLGoldsteinPrice)->ArgsProduct({{10, 20, 50}, {10, 20, 50}});
BENCHMARK(BM_MLSLHartmann3)->ArgsProduct({{10, 20, 50}, {10, 20, 50}});

}  // namespace
}  // namespace numerics
}  // namespace principia
