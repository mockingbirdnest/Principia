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
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "graphics/colours.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/apsides.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/solar_system.hpp"
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
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_apsides;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
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
using namespace principia::testing_utilities::_numerics_matchers;

// A minimum bounding rectangle for a set of values of the eccentricity vector.
struct ExpectedEccentricityVectorRange {
  ApproximateQuantity<double> min_e_cos_ω;
  ApproximateQuantity<double> max_e_cos_ω;
  ApproximateQuantity<double> min_e_sin_ω;
  ApproximateQuantity<double> max_e_sin_ω;
};

struct GeopotentialTruncation {
  // The geopotential truncation used.
  int max_degree;
  bool zonal_only;

  // Expectations.  All values are checked with IsNear.

  // The Euclidean norm of the change in the (e cos ω, e sin ω) vector between
  // the beginning and the end of the first cycle of the repeat ground track
  // orbit.
  ApproximateQuantity<double> first_cycle_eccentricity_vector_drift;

  // An expectation for the periodic repeat ground track behaviour discounting
  // any effects that are periodic over one orbit.
  // Bounds the value of the eccentricity vector at the descending node of every
  // orbit of the first cycle.
  ExpectedEccentricityVectorRange first_cycle_descending_nodes;

  // An expectation for the behaviour over `long_term`.  Only the first
  // descending node of each period is bounded, so that the periodic component
  // tested above is ignored.
  ExpectedEccentricityVectorRange cycle_ends;
  // The duration for the above expectation.
  Time long_term;

  // A string describing the truncation.
  std::string DegreeAndOrder() const {
    return absl::StrCat(max_degree, "×", zonal_only ? 0 : max_degree);
  }
};

std::ostream& operator<<(std::ostream& o, GeopotentialTruncation truncation) {
  return o << truncation.DegreeAndOrder();
}

class LunarOrbitTest : public ::testing::TestWithParam<GeopotentialTruncation> {
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
std::array<GeopotentialTruncation, 7> const geopotential_truncations = {
#else
std::array<GeopotentialTruncation, 5> const geopotential_truncations = {
#endif
    {{
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
         .max_degree = 50,
         .zonal_only = false,
         .first_cycle_eccentricity_vector_drift = 0.00018,
         .first_cycle_descending_nodes = {-0.0055, +0.0051, -0.027, -0.018},
         .cycle_ends = {+0.0026, +0.0037, -0.021, -0.0200},
         .long_term = 10 * JulianYear,
     },
     {
#endif
         .max_degree = 30,
         .zonal_only = true,
         .first_cycle_eccentricity_vector_drift = 0.00032_(1),
         .first_cycle_descending_nodes = {-0.0058_(1), +0.0048_(1), -0.027_(1), -0.018_(1)},
         .cycle_ends = {+0.0019_(1), +0.0050_(1), -0.022_(1), -0.0190_(1)},
         .long_term = 10 * JulianYear,
     },
     {
         .max_degree = 30,
         .zonal_only = false,
         .first_cycle_eccentricity_vector_drift = 0.00032_(1),
         .first_cycle_descending_nodes = {-0.0058_(1), +0.0048_(1), -0.027_(1), -0.018_(1)},
         .cycle_ends = {+0.0019_(1), +0.0050_(1), -0.022_(1), -0.0190_(1)},
         .long_term = 1 * JulianYear,
     },
     {
         .max_degree = 25,
         .zonal_only = false,
         .first_cycle_eccentricity_vector_drift = 0.00110_(1),
         .first_cycle_descending_nodes = {-0.0060_(1), +0.0044_(1), -0.027_(1), -0.018_(1)},
         .cycle_ends = {-0.0017_(1), +0.0089_(1), -0.021_(1), -0.0110_(1)},
         .long_term = 1 * JulianYear,
     },
     {
         .max_degree = 20,
         .zonal_only = false,
         .first_cycle_eccentricity_vector_drift = 0.00130_(1),
         .first_cycle_descending_nodes = {-0.0064_(1), +0.0045_(1), -0.028_(1), -0.018_(1)},
         .cycle_ends = {-0.0030_(1), +0.0100_(1), -0.021_(1), -0.0083_(1)},
         .long_term = 1 * JulianYear,
     },
     {
         .max_degree = 10,
         .zonal_only = false,
         .first_cycle_eccentricity_vector_drift = 0.00370_(1),
         .first_cycle_descending_nodes = {-0.0091_(1), +0.0036_(1), -0.028_(1), -0.018_(1)},
         .cycle_ends = {-0.0160_(1), +0.0210_(1), -0.021_(1), +0.0160_(1)},
         .long_term = 1 * JulianYear,
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
     },
     {
         .max_degree = 50,
         .zonal_only = true,
         .first_cycle_eccentricity_vector_drift = 0.00098,
         .first_cycle_descending_nodes = {+0.0038, +0.0040, -0.022, -0.021},
         .cycle_ends = {-0.0047, +0.0040, -0.025, -0.0170},
         .long_term = 1 * JulianYear,
#endif
     }},
};

INSTANTIATE_TEST_SUITE_P(
    DISABLED_TruncatedSelenopotentials,
    LunarOrbitTest,
    ::testing::ValuesIn(geopotential_truncations));

TEST_P(LunarOrbitTest, NearCircularRepeatGroundTrackOrbit) {
  Time const integration_step = 10 * Second;
  LOG(INFO) << "Using a " << GetParam() << " selenopotential field";

  Logger logger(
      SOLUTION_DIR / "mathematica" /
          absl::StrCat(
              "lunar_orbit_", GetParam().DegreeAndOrder(), ".generated.wl"),
      /*make_unique=*/false);

  // We work with orbit C from [RL06].

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

  logger.Set("tu", TU, ExpressIn(Second));
  logger.Set("lu", LU, ExpressIn(Metre));

  Time const cycle = 2 * π * TU;
  int const revolutions_per_cycle = 328;

  // Initial conditions and elements from table 2 of [RL06].
  Length const x0 = -4.498948742093e-03 * LU;
  Length const y0 = -1.731769313131e-03 * LU;
  Length const z0 =  0 * LU;
  Speed const  u0 = -6.203996010078e-02 * (LU / TU);
  Speed const  v0 =  7.000280770869e-02 * (LU / TU);
  Speed const  w0 =  1.588813067177e+00 * (LU / TU);

  DegreesOfFreedom<LunarSurface> const lunar_initial_state = {
      LunarSurface::origin + Displacement<LunarSurface>({x0, y0, z0}),
      Velocity<LunarSurface>({u0, v0, w0})};

  EXPECT_OK(ephemeris_->Prolong(J2000));
  DegreesOfFreedom<ICRS> const initial_state =
      lunar_frame_.FromThisFrameAtTime(J2000)(lunar_initial_state);

  {
    Length const a0 = +1.861791339407e+03 * Kilo(Metre);
    double const e0 = +2.110475283361e-02;
    Angle const  i0 = +9.298309294740e+01 * Degree;
    Angle const  ω0 = -7.839337618501e+01 * Degree;
    Angle const  Ω0 = -1.589469097527e+02 * Degree;

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
                RelativeErrorFrom(a0, IsNear(9.0e-4_(1))));
    EXPECT_THAT(*initial_osculating.eccentricity,
                RelativeErrorFrom(e0, IsNear(1.4e-10_(1))));
    EXPECT_THAT(initial_osculating.inclination,
                RelativeErrorFrom(i0, IsNear(9.7e-9_(1))));
    EXPECT_THAT(*initial_osculating.argument_of_periapsis,
                RelativeErrorFrom(2 * π * Radian + ω0, IsNear(2.0e-11_(1))));
    EXPECT_THAT(initial_osculating.longitude_of_ascending_node,
                RelativeErrorFrom(2 * π * Radian + Ω0, IsNear(4.7e-13_(1))));
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

  EXPECT_OK(
      ephemeris_->FlowWithFixedStep(J2000 + GetParam().long_term, *instance));

  // To find the nodes, we need to convert the trajectory to a reference frame
  // whose xy plane is the Moon's equator.
  DiscreteTrajectory<Selenocentric> first_cycle_selenocentric_trajectory;
  DiscreteTrajectory<Selenocentric> selenocentric_trajectory;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    if (time - J2000 <= cycle) {
      EXPECT_OK(first_cycle_selenocentric_trajectory.Append(
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

  
  std::vector<EccentricityVector> first_cycle_osculating_eccentricity_vector;
  for (Instant t = J2000; t <= J2000 + cycle; t += cycle / 50'000) {
    auto const elements = KeplerOrbit<Selenocentric>(
            *moon_,
            satellite_,
            selenocentric_trajectory.EvaluateDegreesOfFreedom(t) -
                selenocentre_,
            t).elements_at_epoch();

    Angle const& ω = *elements.argument_of_periapsis;
    double const& e = *elements.eccentricity;
    first_cycle_osculating_eccentricity_vector.emplace_back(e, ω);
  }

  auto const first_cycle_mean_elements = OrbitalElements::ForTrajectory(
      first_cycle_selenocentric_trajectory, *moon_, MasslessBody{});
  ASSERT_OK(first_cycle_mean_elements);
  auto first_cycle_eccentricity_vector_graph =
      first_cycle_mean_elements->PlotEccentricityVector(
          200,
          200,
          /*background=*/Opaque(xkcd::white),
          /*axis_colour=*/xkcd::black,
          /*line_colour=*/xkcd::red);
  first_cycle_eccentricity_vector_graph.ListPointPlot(
      first_cycle_osculating_eccentricity_vector, xkcd::blue);
  EXPECT_GOLDEN_GRAPH(first_cycle_eccentricity_vector_graph,
                      "rl06c_first_cycle");

  DistinguishedPoints<Selenocentric> ascending_nodes;
  DistinguishedPoints<Selenocentric> descending_nodes;
  EXPECT_OK(ComputeNodes(selenocentric_trajectory,
                         selenocentric_trajectory.begin(),
                         selenocentric_trajectory.end(),
                         /*t_max=*/InfiniteFuture,
                         /*north=*/Vector<double, Selenocentric>({0, 0, 1}),
                         /*max_points=*/std::numeric_limits<int>::max(),
                         ascending_nodes,
                         descending_nodes));

  struct Nodes {
    std::string_view const name;
    DistinguishedPoints<Selenocentric> const& points;
  };

  struct Apsides {
    std::string_view const name;
    DistinguishedPoints<ICRS> const& points;
  };

  std::vector<EccentricityVector> descending_node_eccentricity_vector;

  for (auto const& [time, degrees_of_freedom] : descending_nodes) {
    auto const elements = KeplerOrbit<Selenocentric>(
            *moon_,
            satellite_,
            selenocentric_trajectory.EvaluateDegreesOfFreedom(time) -
                selenocentre_,
            time).elements_at_epoch();
      descending_node_eccentricity_vector
          .emplace_back(*elements.eccentricity,*elements.argument_of_periapsis);
  }

  {
    auto const e0 = descending_node_eccentricity_vector[0];
    auto const e1 = descending_node_eccentricity_vector[revolutions_per_cycle];
    EXPECT_THAT(
        Sqrt(Pow<2>(e0.e_cos_ω - e1.e_cos_ω) + Pow<2>(e0.e_sin_ω - e1.e_sin_ω)),
        IsNear(GetParam().first_cycle_eccentricity_vector_drift));
  }

  {
    Interval<double> e_cos_ω;
    Interval<double> e_sin_ω;
    for (auto const& eccentricity_vector :
         descending_node_eccentricity_vector) {
      e_cos_ω.Include(eccentricity_vector.e_cos_ω);
      e_sin_ω.Include(eccentricity_vector.e_sin_ω);
    }
    EXPECT_THAT(e_cos_ω.min,
                IsNear(GetParam().first_cycle_descending_nodes.min_e_cos_ω));
    EXPECT_THAT(e_cos_ω.max,
                IsNear(GetParam().first_cycle_descending_nodes.max_e_cos_ω));
    EXPECT_THAT(e_sin_ω.min,
                IsNear(GetParam().first_cycle_descending_nodes.min_e_sin_ω));
    EXPECT_THAT(e_sin_ω.max,
                IsNear(GetParam().first_cycle_descending_nodes.max_e_sin_ω));
  }

  auto const long_term_mean_elements = OrbitalElements::ForTrajectory(
      selenocentric_trajectory, *moon_, MasslessBody{});
  ASSERT_OK(long_term_mean_elements);
  auto long_term_eccentricity_vector_graph =
      long_term_mean_elements->PlotEccentricityVector(
          200,
          200,
          /*background=*/Opaque(xkcd::white),
          /*axis_colour=*/xkcd::black,
          /*line_colour=*/xkcd::red);
  long_term_eccentricity_vector_graph.ListPointPlot(
      descending_node_eccentricity_vector, xkcd::black);

  {
    Interval<double> e_cos_ω;
    Interval<double> e_sin_ω;
    for (int revolution = 0;
         revolution < descending_node_eccentricity_vector.size();
         revolution += revolutions_per_cycle) {
      auto const eccentricity_vector =
          descending_node_eccentricity_vector[revolution];
      e_cos_ω.Include(eccentricity_vector.e_cos_ω);
      e_sin_ω.Include(eccentricity_vector.e_sin_ω);
      long_term_eccentricity_vector_graph.ListPointPlot(
          std::array{eccentricity_vector}, xkcd::green);
    }
    EXPECT_THAT(e_cos_ω.min,
                IsNear(GetParam().cycle_ends.min_e_cos_ω));
    EXPECT_THAT(e_cos_ω.max,
                IsNear(GetParam().cycle_ends.max_e_cos_ω));
    EXPECT_THAT(e_sin_ω.min,
                IsNear(GetParam().cycle_ends.min_e_sin_ω));
    EXPECT_THAT(e_sin_ω.max,
                IsNear(GetParam().cycle_ends.max_e_sin_ω));
  }
  EXPECT_GOLDEN_GRAPH(long_term_eccentricity_vector_graph,
                      "rl06c_long_term");
}

#endif

}  // namespace astronomy
}  // namespace principia
