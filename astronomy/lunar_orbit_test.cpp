#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "glog/logging.h"
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
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace astronomy {

using ::testing::Lt;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
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
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics;

// A minimum bounding rectangle for a set of values of the eccentricity vector.
struct EccentricityVectorRange {
  double min_e_cos_ω = +Infinity<double>;
  double max_e_cos_ω = -Infinity<double>;
  double min_e_sin_ω = +Infinity<double>;
  double max_e_sin_ω = -Infinity<double>;
};

struct GeopotentialTruncation {
  // The geopotential truncation used.
  int max_degree;
  bool zonal_only;

  // Expectations.  All values are checked with IsNear.

  // The Euclidean norm of the change in the (e cos ω, e sin ω) vector between
  // the beginning and the end of the first period of the repeat ground track
  // orbit.
  double first_period_eccentricity_vector_drift;

  // An expectation for the periodic repeat ground track behaviour discounting
  // any effects that are periodic over one orbit.
  // Bounds the value of the eccentricity vector at the descending node of every
  // orbit of the first period.
  EccentricityVectorRange first_period_descending_nodes;

  // An expectation for the behaviour over periods.  Only the first
  // descending node of each period is bounded, so that the periodic component
  // tested above is ignored.
  EccentricityVectorRange period_ends;
  // The number of periods for the above expectation.
  int periods;

  // A string describing the truncation.
  std::string DegreeAndOrder() const {
    return absl::StrCat(max_degree, "x", zonal_only ? 0 : max_degree);
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
    google::LogToStderr();
  }

  // This Moon-centred, Moon-fixed reference frame has the x axis pointing
  // towards the Earth, and the y axis in the direction of the velocity of the
  // Earth, see figure 1. of [RL06].
  using LunarSurface = Frame<struct LunarSurfaceTag, Arbitrary>;

  // This reference frame is non-rotating, with its origin at the selenocentre.
  // The axes are those of LunarSurface at J2000.
  // Note that this frame is not actually inertial, but we want to use it with
  // `KeplerOrbit`.  Perhaps we should have a concept of non-rotating, and
  // `KeplerOrbit` should check that; this is good enough for a test.
  using Selenocentric = Frame<struct SelenocentricTag, Inertial>;

  // We do not use a `BodyCentredNonRotatingReferenceFrame` since that would use
  // ICRS axes.
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
constexpr std::array<GeopotentialTruncation, 6> geopotential_truncations = {
#else
constexpr std::array<GeopotentialTruncation, 4> geopotential_truncations = {
#endif
    {{
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
         .max_degree = 50,
         .zonal_only = false,
         .first_period_eccentricity_vector_drift = 0.00018,
         .first_period_descending_nodes = {-0.0055, +0.0051, -0.027, -0.018},
         .period_ends = {+0.0026, +0.0037, -0.021, -0.0200},
         .periods = 10,
     },
     {
#endif
         .max_degree = 30,
         .zonal_only = false,
         .first_period_eccentricity_vector_drift = 0.00032,
         .first_period_descending_nodes = {-0.0058, +0.0048, -0.027, -0.018},
         .period_ends = {+0.0019, +0.0050, -0.022, -0.0190},
         .periods = 28,
     },
     {
         .max_degree = 25,
         .zonal_only = false,
         .first_period_eccentricity_vector_drift = 0.00110,
         .first_period_descending_nodes = {-0.0060, +0.0044, -0.027, -0.018},
         .period_ends = {-0.0017, +0.0089, -0.021, -0.0110},
         .periods = 28,
     },
     {
         .max_degree = 20,
         .zonal_only = false,
         .first_period_eccentricity_vector_drift = 0.00130,
         .first_period_descending_nodes = {-0.0064, +0.0045, -0.028, -0.018},
         .period_ends = {-0.0030, +0.0100, -0.021, -0.0083},
         .periods = 28,
     },
     {
         .max_degree = 10,
         .zonal_only = false,
         .first_period_eccentricity_vector_drift = 0.00370,
         .first_period_descending_nodes = {-0.0091, +0.0036, -0.028, -0.018},
         .period_ends = {-0.0160, +0.0210, -0.021, +0.0160},
         .periods = 28,
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
     },
     {
         .max_degree = 50,
         .zonal_only = true,
         .first_period_eccentricity_vector_drift = 0.00098,
         .first_period_descending_nodes = {+0.0038, +0.0040, -0.022, -0.021},
         .period_ends = {-0.0047, +0.0040, -0.025, -0.0170},
         .periods = 28,
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
  EXPECT_THAT(RelativeError(TU, TU_rl), IsNear(1.4e-3_(1)));
  EXPECT_THAT(RelativeError(LU, LU_rl), IsNear(9.0e-4_(1)));

  logger.Set("tu", TU, ExpressIn(Second));
  logger.Set("lu", LU, ExpressIn(Metre));

  Time const period = 2 * π * TU;
  int const orbits_per_period = 328;

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
    EXPECT_THAT(RelativeError(*initial_osculating.semimajor_axis, a0),
                IsNear(9.0e-4_(1)));
    EXPECT_THAT(RelativeError(*initial_osculating.eccentricity, e0),
                IsNear(1.4e-10_(1)));
    EXPECT_THAT(RelativeError(initial_osculating.inclination, i0),
                IsNear(9.7e-9_(1)));
    EXPECT_THAT(RelativeError(*initial_osculating.argument_of_periapsis,
                              2 * π * Radian + ω0),
                IsNear(2.0e-11_(1)));
    EXPECT_THAT(RelativeError(initial_osculating.longitude_of_ascending_node,
                              2 * π * Radian + Ω0),
                IsNear(4.7e-13_(1)));
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

  EXPECT_OK(ephemeris_->FlowWithFixedStep(J2000 + GetParam().periods * period,
                                          *instance));

  // To find the nodes, we need to convert the trajectory to a reference frame
  // whose xy plane is the Moon's equator.
  DiscreteTrajectory<LunarSurface> surface_trajectory;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    EXPECT_OK(surface_trajectory.Append(
        time, lunar_frame_.ToThisFrameAtTime(time)(degrees_of_freedom)));
  }

  for (Instant t = J2000; t <= J2000 + 2 * period; t += period / 50'000) {
    auto const elements = KeplerOrbit<Selenocentric>(
        *moon_,
        satellite_,
        ToSelenocentric(t)(
            trajectory.EvaluateDegreesOfFreedom(t)) - selenocentre_,
        t).elements_at_epoch();

    logger.Append("times",
                  t - J2000,
                  ExpressIn(Second));
    logger.Append("semimajorAxes",
                  *elements.semimajor_axis,
                  ExpressIn(Metre));
    logger.Append("inclinations",
                  elements.inclination,
                  ExpressIn(Radian));
    logger.Append("eccentricities",
                  *elements.eccentricity);
    logger.Append("arguments",
                  *elements.argument_of_periapsis,
                  ExpressIn(Radian));
    logger.Append("longitudesOfAscendingNodes",
                  elements.longitude_of_ascending_node,
                  ExpressIn(Radian));
    logger.Append("displacements",
                  surface_trajectory.EvaluatePosition(t) - LunarSurface::origin,
                  ExpressIn(Metre));
  }

  DiscreteTrajectory<LunarSurface> ascending_nodes;
  DiscreteTrajectory<LunarSurface> descending_nodes;
  EXPECT_OK(ComputeNodes(surface_trajectory,
                         surface_trajectory.begin(),
                         surface_trajectory.end(),
                         /*t_max=*/InfiniteFuture,
                         /*north=*/Vector<double, LunarSurface>({0, 0, 1}),
                         /*max_points=*/std::numeric_limits<int>::max(),
                         ascending_nodes,
                         descending_nodes));

  DiscreteTrajectory<ICRS> apoapsides;
  DiscreteTrajectory<ICRS> periapsides;
  ComputeApsides(*ephemeris_->trajectory(moon_),
                 trajectory,
                 trajectory.begin(),
                 trajectory.end(),
                 /*t_max=*/InfiniteFuture,
                 /*max_points=*/std::numeric_limits<int>::max(),
                 apoapsides,
                 periapsides);

  struct Nodes {
    std::string_view const name;
    DiscreteTrajectory<LunarSurface> const& trajectory;
  };

  struct Apsides {
    std::string_view const name;
    DiscreteTrajectory<ICRS> const& trajectory;
  };

  std::vector<double> descending_node_eccentricities;
  std::vector<Angle> descending_node_arguments;

  for (auto const& nodes : {Nodes{"ascending", ascending_nodes},
                            Nodes{"descending", descending_nodes}}) {
    for (auto const& [time, degrees_of_freedom] : nodes.trajectory) {
      auto const elements = KeplerOrbit<Selenocentric>(
          *moon_,
          satellite_,
          ToSelenocentric(time)(
              trajectory.EvaluateDegreesOfFreedom(time)) - selenocentre_,
          time).elements_at_epoch();

      logger.Append(absl::StrCat(nodes.name, "NodeTimes"),
                    time - J2000,
                    ExpressIn(Second));
      logger.Append(absl::StrCat(nodes.name, "NodeDisplacements"),
                    degrees_of_freedom.position() - LunarSurface::origin,
                    ExpressIn(Metre));
      logger.Append(absl::StrCat(nodes.name, "NodeArguments"),
                    *elements.argument_of_periapsis,
                    ExpressIn(Radian));
      logger.Append(absl::StrCat(nodes.name, "NodeEccentricities"),
                    *elements.eccentricity);

      if (nodes.name == "descending") {
        descending_node_eccentricities.push_back(*elements.eccentricity);
        descending_node_arguments.push_back(*elements.argument_of_periapsis);
      }
    }
  }

  for (auto const& apsides : {Apsides{"apoapsis", apoapsides},
                              Apsides{"periapsis", periapsides}}) {
    for (auto const& [time, degrees_of_freedom] : apsides.trajectory) {
      logger.Append(absl::StrCat(apsides.name, "Times"),
                    time - J2000,
                    ExpressIn(Second));
      logger.Append(absl::StrCat(apsides.name, "Displacements"),
                    lunar_frame_.ToThisFrameAtTime(time).rigid_transformation()(
                        degrees_of_freedom.position()) -
                        LunarSurface::origin,
                    ExpressIn(Metre));
    }
  }

  {
    auto const e0 = descending_node_eccentricities[0];
    auto const e1 = descending_node_eccentricities[orbits_per_period];
    auto const ω0 = descending_node_arguments[0];
    auto const ω1 = descending_node_arguments[orbits_per_period];
    EXPECT_THAT(RelativeError(GetParam().first_period_eccentricity_vector_drift,
                              Sqrt(Pow<2>(e1 * Cos(ω1) - e0 * Cos(ω0)) +
                                   Pow<2>(e1 * Sin(ω1) - e0 * Sin(ω0)))),
                Lt(0.035));
  }

  {
    EccentricityVectorRange actual_first_period_descending_nodes;
    for (int orbit = 0; orbit < orbits_per_period; ++orbit) {
      auto& actual = actual_first_period_descending_nodes;
      auto const e = descending_node_eccentricities[orbit];
      auto const ω = descending_node_arguments[orbit];
      actual.min_e_cos_ω = std::min(actual.min_e_cos_ω, e * Cos(ω));
      actual.max_e_cos_ω = std::max(actual.max_e_cos_ω, e * Cos(ω));
      actual.min_e_sin_ω = std::min(actual.min_e_sin_ω, e * Sin(ω));
      actual.max_e_sin_ω = std::max(actual.max_e_sin_ω, e * Sin(ω));
    }
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.min_e_cos_ω,
                      actual_first_period_descending_nodes.min_e_cos_ω),
        Lt(0.007));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.max_e_cos_ω,
                      actual_first_period_descending_nodes.max_e_cos_ω),
        Lt(0.012));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.min_e_sin_ω,
                      actual_first_period_descending_nodes.min_e_sin_ω),
        Lt(0.017));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.max_e_sin_ω,
                      actual_first_period_descending_nodes.max_e_sin_ω),
        Lt(0.017));
  }

  {
    EccentricityVectorRange actual_period_ends;
    for (int orbit = 0;
         orbit < descending_nodes.size();
         orbit += orbits_per_period) {
      auto& actual = actual_period_ends;
      auto const e = descending_node_eccentricities[orbit];
      auto const ω = descending_node_arguments[orbit];
      actual.min_e_cos_ω = std::min(actual.min_e_cos_ω, e * Cos(ω));
      actual.max_e_cos_ω = std::max(actual.max_e_cos_ω, e * Cos(ω));
      actual.min_e_sin_ω = std::min(actual.min_e_sin_ω, e * Sin(ω));
      actual.max_e_sin_ω = std::max(actual.max_e_sin_ω, e * Sin(ω));
    }
    EXPECT_THAT(RelativeError(GetParam().period_ends.min_e_cos_ω,
                              actual_period_ends.min_e_cos_ω),
                Lt(0.015));
    EXPECT_THAT(RelativeError(GetParam().period_ends.max_e_cos_ω,
                              actual_period_ends.max_e_cos_ω),
                Lt(0.019));
    EXPECT_THAT(RelativeError(GetParam().period_ends.min_e_sin_ω,
                              actual_period_ends.min_e_sin_ω),
                Lt(0.017));
    EXPECT_THAT(RelativeError(GetParam().period_ends.max_e_sin_ω,
                              actual_period_ends.max_e_sin_ω),
                Lt(0.027));
  }
}

#endif

}  // namespace astronomy
}  // namespace principia
