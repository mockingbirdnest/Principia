#include <algorithm>
#include <array>
#include <filesystem>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

#include "absl/log/check.h"
#include "absl/log/globals.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/orbital_elements.hpp"
#include "astronomy/лидов.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "graphics/colours.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "numerics/angle_reduction.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/apsides.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_view.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/golden_graphs.hpp"  // 🧙 For EXPECT_GOLDEN_GRAPH.
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace astronomy {

using ::testing::Lt;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_лидов;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::graphics::_colours;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::numerics::_angle_reduction;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_apsides;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_discrete_trajectory_view;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_oblate_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

// A minimum bounding rectangle for a set of values of the eccentricity vector.
struct ExpectedEccentricityVectorRange {
  ApproximateQuantity<double> min_e_cos_ω;
  ApproximateQuantity<double> max_e_cos_ω;
  ApproximateQuantity<double> min_e_sin_ω;
  ApproximateQuantity<double> max_e_sin_ω;
};

struct RL06Orbit {
  // Initial conditions and elements from table 2 of [RL06].
  int const revolutions_per_month;
  Length const x₀;
  Length const y₀;
  Length const z₀ = 0 * Metre;
  Speed const  u₀;
  Speed const v₀;
  Speed const w₀;
  Length const a₀;
  double const e₀;
  Angle const i₀;
  Angle const ω₀;
  Angle const Ω₀;
  // Expected relative errors on the initial osculating elements.
  // The error of a₀ relatively large (9.0e-4) and is independent of the orbit:
  // it is the error between our LU and the one from [RL06].  The errors of the
  // eccentricity and angles should be tiny, reflecting differences in
  // evaluation of the geometric computations; they serve as a check on data
  // entry.
  ApproximateQuantity<double> const e₀_error;
  ApproximateQuantity<double> const i₀_error;
  ApproximateQuantity<double> const ω₀_error;
  ApproximateQuantity<double> const Ω₀_error;

  static RL06Orbit A(Length const& LU, Time const& TU) {
    return {
        .revolutions_per_month = 73,
        .x₀ = -1.311519120505e-02  * LU,
        .y₀ = 5.435394815081e-04 * LU,
        .u₀ = -1.281711107594e-02 * (LU / TU),
        .v₀ = -3.056086111584e-01 * (LU / TU),
        .w₀ = 9.077797920947e-01 * (LU / TU),
        .a₀ = 5.046738218681e+03 * Kilo(Metre),
        .e₀ = 2.425326521133e-04,
        .i₀ = 7.063797094157e+01 * Degree,
        .ω₀ = -4.048674547795e+01 * Degree,
        .Ω₀ = 1.776268201967e+02 * Degree,
        .e₀_error = 4.0e-9_(1),
        .i₀_error = 9.9e-15_(1),
        .ω₀_error = 1.2e-10_(1),
        .Ω₀_error = 2.5e-13_(1),
    };
  }

  static RL06Orbit B(Length const& LU, Time const& TU) {
    return {
        .revolutions_per_month = 73,
        .x₀ = -7.660645625403e-03  * LU,
        .y₀ = 5.028162133106e-03 * LU,
        .u₀ = 1.327777508469e-01 * (LU / TU),
        .v₀ = -9.233621714500e-01 * (LU / TU),
        .w₀ = 9.132878998189e-01 * (LU / TU),
        .a₀ = 4.996647749602e+03 * Kilo(Metre),
        .e₀ = 5.384086098625e-01,
        .i₀ = 5.220698531621e+01 * Degree,
        .ω₀ = 8.922084663298e+01 * Degree,
        .Ω₀ = 1.467205983429e+02 * Degree,
        .e₀_error = 1.5e-12_(1),
        .i₀_error = 3.4e-13_(1),
        .ω₀_error = 1.3e-12_(1),
        .Ω₀_error = 3.6e-13_(1),
    };
  }

  static RL06Orbit C(Length const& LU, Time const& TU) {
    return {
        .revolutions_per_month = 328,
        .x₀ = -4.498948742093e-03 * LU,
        .y₀ = -1.731769313131e-03 * LU,
        .u₀ = -6.203996010078e-02 * (LU / TU),
        .v₀ = 7.000280770869e-02 * (LU / TU),
        .w₀ = 1.588813067177e+00 * (LU / TU),
        .a₀ = +1.861791339407e+03 * Kilo(Metre),
        .e₀ = +2.110475283361e-02,
        .i₀ = +9.298309294740e+01 * Degree,
        .ω₀ = -7.839337618501e+01 * Degree,
        .Ω₀ = -1.589469097527e+02 * Degree,
        .e₀_error = 1.4e-10_(1),
        .i₀_error = 9.7e-9_(1),
        .ω₀_error = 2.0e-11_(1),
        .Ω₀_error = 4.7e-13_(1),
    };
  }
};


struct OrbitAndGeopotentialTruncation {
  // The orbit from [RL06], as a function of the lunar length and time units;
  RL06Orbit (&orbit)(Length const& LU, Time const& TU);
  // The geopotential truncation used.
  int max_degree;
  bool zonal_only;

  // Expectations.  All values are checked with IsNear.

  // The Euclidean norm of the change in the (e cos ω, e sin ω) vector between
  // the beginning and the end of the first period of the ground track cycle
  // orbit (which lasts one month).
  ApproximateQuantity<double> first_month_eccentricity_vector_drift;

  // An expectation for the periodic repeat ground track behaviour discounting
  // any effects that are periodic over one orbit.
  // Bounds the value of the eccentricity vector at the descending node of every
  // orbit of the first month.
  ExpectedEccentricityVectorRange first_month_descending_nodes;

  // An expectation for the behaviour over `long_term`.  Only the first
  // descending node of each period is bounded, so that the periodic component
  // tested above is ignored.
  ExpectedEccentricityVectorRange month_ends;
  // The duration for the above expectation.
  Time long_term;

  // The expected status of the integration over `long_term`.
  absl::StatusCode long_term_status = absl::StatusCode::kOk;

  // Override the range of the eccentricity vector plot for the first month.
  // This is used for orbit A, where the osculating eccentricity vector varies
  // much more than the mean one.
  std::optional<double> first_month_e_cos_ω_plot_range;

  std::string_view OrbitName() const {
    if (&orbit == &RL06Orbit::A) {
      return "A";
    } else if (&orbit == &RL06Orbit::B) {
      return "B";
    } else if (&orbit == &RL06Orbit::C) {
      return "C";
    } else {
      LOG(FATAL) << "Unknown orbit " << orbit;
    }
  }

  // A string describing the truncation.
  std::string DegreeAndOrder() const {
    return absl::StrCat(max_degree, "×", zonal_only ? 0 : max_degree);
  }
};

std::ostream& operator<<(std::ostream& o,
                         OrbitAndGeopotentialTruncation param) {
  return o << param.OrbitName() << "_" << param.DegreeAndOrder();
}

class LunarOrbitTest
    : public ::testing::TestWithParam<OrbitAndGeopotentialTruncation> {
 protected:
  LunarOrbitTest()
      : solar_system_2000_([]() {
          SolarSystem<ICRS> result(
              SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
              SOLUTION_DIR / "astronomy" /
                  "sol_initial_state_jd_2451545_000000000.proto.txt");
          result.LimitOblatenessToDegree("Moon", GetParam().max_degree);
          if (GetParam().zonal_only) {
            result.LimitOblatenessToZonal("Moon");
          }
          return result;
        }()),
        ephemeris_(solar_system_2000_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<
                    QuinlanTremaine1990Order12,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*step=*/10 * Minute))),
        moon_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_2000_.massive_body(*ephemeris_, "Moon"))),
        lunar_frame_(ephemeris_.get(), moon_),
        selenocentre_(Selenocentric::origin, Selenocentric::unmoving) {
    absl::SetStderrThreshold(absl::LogSeverityAtLeast::kInfo);
  }

  // This Moon-centred, Moon-fixed reference frame has the x axis pointing
  // towards the Earth, and the y axis in the direction of the velocity of the
  // Earth, see figure 1. of [RL06].
  using LunarSurface = Frame<struct LunarSurfaceTag, Arbitrary>;

  // This reference frame is non-rotating, with its origin at the selenocentre.
  // The axes are those of LunarSurface at J2000.
  using Selenocentric = Frame<struct SelenocentricTag, NonRotating>;

  // We do not use a `BodyCentredNonRotatingReferenceFrame` since that would not
  // have its x axis pointing towards the Earth at J2000.
  RigidMotion<ICRS, Selenocentric> ToSelenocentric(Instant const& t) {
    return RigidMotion<ICRS, Selenocentric>(
        RigidTransformation<ICRS, Selenocentric>(
            ephemeris_->trajectory(moon_)->EvaluatePosition(t),
            Selenocentric::origin,
            OrthogonalMap<LunarSurface, Selenocentric>::Identity() *
                lunar_frame_.ToThisFrameAtTime(J2000).orthogonal_map()),
        /*angular_velocity_of_to_frame=*/ICRS::nonrotating,
        /*velocity_of_to_frame_origin=*/
        ephemeris_->trajectory(moon_)->EvaluateVelocity(t));
  }

  SolarSystem<ICRS> const solar_system_2000_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<OblateBody<ICRS> const*> const moon_;

  BodySurfaceReferenceFrame<ICRS, LunarSurface> const lunar_frame_;
  DegreesOfFreedom<Selenocentric> selenocentre_;

  MasslessBody const satellite_;
};

#if !defined(_DEBUG)

#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
std::array<OrbitAndGeopotentialTruncation, 8> const geopotential_truncations = {
#else
std::array<OrbitAndGeopotentialTruncation, 7> const geopotential_truncations = {
#endif
    {{
         .orbit = RL06Orbit::A,
         .max_degree = 30,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.000'061_(1),
         .first_month_descending_nodes =
             {-0.51_(1), +0.000'043_(1), -0.51_(1), -0.000'062_(1)},
         .month_ends = {-0.50_(1), -0.000'50_(1), -0.50_(1), -0.000'33_(1)},
         // We find a collision somewhere in March 2004; [RL06] has it
         // after 5.39 years.
         .long_term = 4.25 * JulianYear,
         .long_term_status = absl::StatusCode::kOutOfRange,
         .first_month_e_cos_ω_plot_range = 14e-4,
     },
     {
         .orbit = RL06Orbit::B,
         .max_degree = 30,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.0098_(1),
         .first_month_descending_nodes =
             {-0.091_(1), +0.10_(1), +0.42_(1), +0.56_(1)},
         .month_ends = {-0.090_(1), +0.096_(1), +0.42_(1), +0.56_(1)},
         .long_term = 10 * JulianYear,
     },
     {
// This test was used to determine that 30 is an appropriate maximum degree for
// the geopotential, by comparing the appearance of the graphs for orbit C (the
// lowest of the three) in the 50×50 and 𝑛×𝑛 selenopotentials.
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
         .orbit = RL06Orbit::C,
         .max_degree = 50,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.00018_(1),
         .first_month_descending_nodes =
             {-0.0071_(1), +0.0067_(1), -0.029_(1), -0.016_(1)},
         .month_ends = {+0.0020_(1), +0.0052_(1), -0.022_(1), -0.018_(1)},
         .long_term = 10 * JulianYear,
     },
     {
#endif
         .orbit = RL06Orbit::C,
         .max_degree = 30,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.00032_(1),
         .first_month_descending_nodes =
             {-0.0080_(1), +0.0074_(1), -0.029_(1), -0.014_(1)},
         .month_ends = {+0.00062_(1), +0.0055_(1), -0.022_(1), -0.018_(1)},
         .long_term = 10 * JulianYear,
     },
     {
         .orbit = RL06Orbit::C,
         .max_degree = 25,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.0011_(1),
         .first_month_descending_nodes =
             {-0.011_(1), +0.0044_(1), -0.027_(1), -0.0093_(1)},
         .month_ends = {-0.0017_(1), +0.0037_(1), -0.021_(1), -0.011_(1)},
         .long_term = 1 * JulianYear,
     },
     {
         .orbit = RL06Orbit::C,
         .max_degree = 20,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.0013_(1),
         .first_month_descending_nodes =
             {-0.012_(1), +0.0045_(1), -0.028_(1), -0.0069_(1)},
         .month_ends = {-0.0030_(1), +0.0036_(1), -0.021_(1), -0.0088_(1)},
         .long_term = 1 * JulianYear,
     },
     {
         .orbit = RL06Orbit::C,
         .max_degree = 10,
         .zonal_only = false,
         .first_month_eccentricity_vector_drift = 0.0037_(1),
         .first_month_descending_nodes =
             {-0.027_(1), +0.0036_(1), -0.028_(1), +0.014_(1)},
         .month_ends = {-0.016_(1), +0.0036_(1), -0.021_(1), +0.013_(1)},
         .long_term = 1 * JulianYear,
     },
     // Figure 13 from [RL06] compares long-term evolutions in the 50×0 and
     // 50×50 field for an orbit in the same family as orbit C.  The paper does
     // not give that orbit (for which the 50×0 field results in a collision),
     // but we do a similar comparison with Orbit C itself in the 30×0 and 30×30
     // fields.
     {
         .orbit = RL06Orbit::C,
         .max_degree = 30,
         .zonal_only = true,
         .first_month_eccentricity_vector_drift = 0.0011_(1),
         .first_month_descending_nodes =
             {-0.0053_(1), +0.0050_(1), -0.025_(1), -0.013_(1)},
         .month_ends = {-0.0052_(1), +0.0050_(1), -0.025_(1), -0.013_(1)},
         .long_term = 10 * JulianYear,
     }},
};

INSTANTIATE_TEST_SUITE_P(
    DISABLED_TruncatedSelenopotentials,
    LunarOrbitTest,
    ::testing::ValuesIn(geopotential_truncations));

TEST_P(LunarOrbitTest, OrbitalElements) {
  Time const integration_step = 10 * Second;
  LOG(INFO) << "Orbit " << GetParam().OrbitName() << " using a "
            << GetParam().DegreeAndOrder() << " selenopotential field";

  // The length and time units LU and TU are such that, in an idealized
  // Earth-Moon system, the Earth-Moon distance is 1 LU and the angular
  // frequency of the body-fixed moon frame is θ′ = 1 rad / TU.

  // In order to best reproduce the results of the paper, we choose our TU such
  // that the rotational period of the moon is TU, and our LU such that the
  // Moon's gravitational parameter has the same value in LU³/TU² as the
  // gravitational parameter used in the paper.
  // With cartesian initial conditions in the surface frame, these two
  // properties ensure that the initial osculating lunar orbit has the same
  // orientation, eccentricity, and anomaly.

  // The _rl values are the ones from table 1 of [RL06].
  Length const LU_rl = 384'400 * Kilo(Metre);
  Time const TU_rl = 375'190.258663027 * Second;
  GravitationalParameter const GM_rl =
      4'902.801076 * (Pow<3>(Kilo(Metre)) / Pow<2>(Second));

  Time const TU = 1 * Radian / moon_->angular_velocity().Norm();
  Length const LU = Cbrt((moon_->gravitational_parameter() * Pow<2>(TU)) /
                         (GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl))));
  EXPECT_THAT(moon_->gravitational_parameter() / (Pow<3>(LU) / Pow<2>(TU)),
              AlmostEquals(GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl)), 1));
  EXPECT_THAT(TU, RelativeErrorFrom(TU_rl, IsNear(1.4e-3_(1))));
  EXPECT_THAT(LU, RelativeErrorFrom(LU_rl, IsNear(9.0e-4_(1))));

  Time const month = 2 * π * TU;

  auto const orbit = GetParam().orbit(LU, TU);

  DegreesOfFreedom<LunarSurface> const lunar_initial_state = {
      LunarSurface::origin +
          Displacement<LunarSurface>({orbit.x₀, orbit.y₀, orbit.z₀}),
      Velocity<LunarSurface>({orbit.u₀, orbit.v₀, orbit.w₀})};

  EXPECT_OK(ephemeris_->Prolong(J2000));
  DegreesOfFreedom<ICRS> const initial_state =
      lunar_frame_.FromThisFrameAtTime(J2000)(lunar_initial_state);

  {
    KeplerianElements<Selenocentric> const initial_osculating =
        KeplerOrbit<Selenocentric>(
            *moon_,
            satellite_,
            ToSelenocentric(J2000)(initial_state) - selenocentre_,
            J2000).elements_at_epoch();
    // The relative error on the semimajor axis is the same as the relative
    // error on our LU with respect to the one in the paper: the semimajor axis
    // has the same value in LU.
    EXPECT_THAT(*initial_osculating.semimajor_axis,
                RelativeErrorFrom(orbit.a₀, IsNear(9.0e-4_(1))));
    EXPECT_THAT(*initial_osculating.eccentricity,
                RelativeErrorFrom(orbit.e₀, IsNear(orbit.e₀_error)));
    EXPECT_THAT(initial_osculating.inclination,
                RelativeErrorFrom(orbit.i₀, IsNear(orbit.i₀_error)));
    EXPECT_THAT(
        *initial_osculating.argument_of_periapsis,
                RelativeErrorFrom(ReduceAngle<0.0, 2 * π>(orbit.ω₀),
                                  IsNear(orbit.ω₀_error)));
    EXPECT_THAT(
        initial_osculating.longitude_of_ascending_node,
                RelativeErrorFrom(ReduceAngle<0.0, 2 * π>(orbit.Ω₀),
                                  IsNear(orbit.Ω₀_error)));
  }

  DiscreteTrajectory<ICRS> trajectory;
  EXPECT_OK(trajectory.Append(J2000, initial_state));
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          integration_step));
  EXPECT_THAT(
      ephemeris_->FlowWithFixedStep(J2000 + GetParam().long_term, *instance),
      StatusIs(GetParam().long_term_status));

  // To find the nodes, we need to convert the trajectory to a reference frame
  // whose xy plane is the Moon's equator.
  DiscreteTrajectory<Selenocentric> first_month_selenocentric_trajectory;
  DiscreteTrajectory<Selenocentric> selenocentric_trajectory;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    if (time <= J2000 + month) {
      EXPECT_OK(first_month_selenocentric_trajectory.Append(
          time, ToSelenocentric(time)(degrees_of_freedom)))
    }
    EXPECT_OK(selenocentric_trajectory.Append(
        time, ToSelenocentric(time)(degrees_of_freedom)))
  }

  struct EccentricityVector {
    EccentricityVector(double const& e, Angle const& ω) {
      auto const [sin_ω, cos_ω] = SinCos(ω);
      e_sin_ω = e * sin_ω;
      e_cos_ω = e * cos_ω;
    }
    double e_cos_ω;
    double e_sin_ω;
  };

  std::vector<EccentricityVector> first_month_osculating_eccentricity_vector;
  for (Instant t = J2000; t <= J2000 + month; t += month / 50'000) {
    auto const elements = KeplerOrbit<Selenocentric>(
            *moon_,
            satellite_,
            selenocentric_trajectory.EvaluateDegreesOfFreedom(t) -
                selenocentre_,
            t).elements_at_epoch();

    Angle const& ω = *elements.argument_of_periapsis;
    double const& e = *elements.eccentricity;
    first_month_osculating_eccentricity_vector.emplace_back(e, ω);
  }

  auto const first_month_mean_elements = OrbitalElements::ForTrajectory(
      first_month_selenocentric_trajectory, *moon_, MasslessBody{});
  ASSERT_OK(first_month_mean_elements);
  auto first_month_eccentricity_vector_graph =
      first_month_mean_elements->PlotEccentricityVector(
          /*width=*/200,
          /*height=*/200,
          /*background=*/Opaque(xkcd::white),
          /*axis_colour=*/xkcd::black,
          /*line_colour=*/xkcd::red,
          /*e_cos_ω_measure=*/GetParam().first_month_e_cos_ω_plot_range);
  first_month_eccentricity_vector_graph.ListPointPlot(
      first_month_osculating_eccentricity_vector, xkcd::blue);
  // This corresponds to the left-hand side of [RL06] figures 9–11 for orbits
  // A–C, with the addition of the mean elements in red.
  EXPECT_GOLDEN_GRAPH(first_month_eccentricity_vector_graph,
                      "first_month");

  DistinguishedPoints<Selenocentric> ascending_nodes;
  DistinguishedPoints<Selenocentric> descending_nodes;
  EXPECT_OK(ComputeNodes(DiscreteTrajectoryView(&selenocentric_trajectory),
                         /*north=*/Vector<double, Selenocentric>({0, 0, 1}),
                         /*max_points=*/std::numeric_limits<int>::max(),
                         ascending_nodes,
                         descending_nodes));
  struct Nodes {
    std::string_view const name;
    DistinguishedPoints<Selenocentric> const& points;
  };

  std::vector<EccentricityVector> descending_node_eccentricity_vector;

  for (auto const& [time, degrees_of_freedom] : descending_nodes) {
    auto const elements = KeplerOrbit<Selenocentric>(
            *moon_,
            satellite_,
            selenocentric_trajectory.EvaluateDegreesOfFreedom(time) -
                selenocentre_,
            time).elements_at_epoch();
    descending_node_eccentricity_vector.emplace_back(
        *elements.eccentricity, *elements.argument_of_periapsis);
  }

  {
    auto const e0 = descending_node_eccentricity_vector[0];
    auto const e1 =
        descending_node_eccentricity_vector[orbit.revolutions_per_month];
    EXPECT_THAT(
        Sqrt(Pow<2>(e0.e_cos_ω - e1.e_cos_ω) + Pow<2>(e0.e_sin_ω - e1.e_sin_ω)),
        IsNear(GetParam().first_month_eccentricity_vector_drift));
  }

  {
    Interval<double> first_month_e_cos_ω;
    Interval<double> first_month_e_sin_ω;
    for (auto const& eccentricity_vector :
         descending_node_eccentricity_vector) {
      first_month_e_cos_ω.Include(eccentricity_vector.e_cos_ω);
      first_month_e_sin_ω.Include(eccentricity_vector.e_sin_ω);
    }
    EXPECT_THAT(first_month_e_cos_ω.min,
                IsNear(GetParam().first_month_descending_nodes.min_e_cos_ω));
    EXPECT_THAT(first_month_e_cos_ω.max,
                IsNear(GetParam().first_month_descending_nodes.max_e_cos_ω));
    EXPECT_THAT(first_month_e_sin_ω.min,
                IsNear(GetParam().first_month_descending_nodes.min_e_sin_ω));
    EXPECT_THAT(first_month_e_sin_ω.max,
                IsNear(GetParam().first_month_descending_nodes.max_e_sin_ω));
  }

  auto const long_term_mean_elements = OrbitalElements::ForTrajectory(
      selenocentric_trajectory, *moon_, MasslessBody{});
  ASSERT_OK(long_term_mean_elements);
  auto long_term_eccentricity_vector_graph =
      long_term_mean_elements->PlotEccentricityVector(
          /*width=*/200,
          /*height=*/200,
          /*background=*/Opaque(xkcd::white),
          /*axis_colour=*/xkcd::black,
          /*line_colour=*/xkcd::red);
  long_term_eccentricity_vector_graph.ListPointPlot(
      descending_node_eccentricity_vector, xkcd::black);

  {
    Interval<double> long_term_e_cos_ω;
    Interval<double> long_term_e_sin_ω;
    for (int revolution = 0;
         revolution < descending_node_eccentricity_vector.size();
         revolution += orbit.revolutions_per_month) {
      auto const eccentricity_vector =
          descending_node_eccentricity_vector[revolution];
      long_term_e_cos_ω.Include(eccentricity_vector.e_cos_ω);
      long_term_e_sin_ω.Include(eccentricity_vector.e_sin_ω);
      long_term_eccentricity_vector_graph.ListPointPlot(
          std::array{eccentricity_vector}, xkcd::green);
    }
    EXPECT_THAT(long_term_e_cos_ω.min,
                IsNear(GetParam().month_ends.min_e_cos_ω));
    EXPECT_THAT(long_term_e_cos_ω.max,
                IsNear(GetParam().month_ends.max_e_cos_ω));
    EXPECT_THAT(long_term_e_sin_ω.min,
                IsNear(GetParam().month_ends.min_e_sin_ω));
    EXPECT_THAT(long_term_e_sin_ω.max,
                IsNear(GetParam().month_ends.max_e_sin_ω));
  }
  // This corresponds to the right-hand side of [RL06] figures 9–11 for orbits
  // A–C, with the addition of the mean elements in red, and with the points at
  // the end of each ground track cycle shown in green instead of black.
  EXPECT_GOLDEN_GRAPH(long_term_eccentricity_vector_graph,
                      "long_term");

  EXPECT_GOLDEN_GRAPH(ЛидовGraph(*long_term_mean_elements,
                                 /*width=*/200,
                                 /*height=*/200,
                                 /*background=*/Opaque(xkcd::black),
                                 /*region_boundary_colour=*/xkcd::white,
                                 /*inclination_colour=*/xkcd::lavender,
                                 /*eccentricity_colour=*/xkcd::cornflower,
                                 /*лидов_parameter_colour=*/xkcd::rose_red,
                                 ЛидовGrid::MaxEccentricityMinInclination),
                      "лидов");
}

#endif

}  // namespace astronomy
}  // namespace principia
