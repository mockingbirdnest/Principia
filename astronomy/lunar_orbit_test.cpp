
#include <algorithm>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "base/file.hpp"
#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/apsides.hpp"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using astronomy::ICRS;
using astronomy::J2000;
using base::dynamic_cast_not_null;
using base::not_null;
using base::OFStream;
using geometry::AngularVelocity;
using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::ComputeApsides;
using physics::ComputeNodes;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MasslessBody;
using physics::OblateBody;
using physics::RelativeDegreesOfFreedom;
using physics::RigidMotion;
using physics::RigidTransformation;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::ArcSin;
using quantities::Cbrt;
using quantities::Cos;
using quantities::Sin;
using quantities::GravitationalParameter;
using quantities::Infinity;
using quantities::Length;
using quantities::Pow;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::ArcMinute;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::Slope;
using testing_utilities::operator""_⑴;

namespace astronomy {

// A minimum bounding rectangle for a set of values of the eccentricity vector.
struct EccentricityVectorRange {
  double min_e_cos_ω = +Infinity<double>();
  double max_e_cos_ω = -Infinity<double>();
  double min_e_sin_ω = +Infinity<double>();
  double max_e_sin_ω = -Infinity<double>();
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
      : solar_system_2000_([this]() {
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
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/10 * Minute))),
        moon_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_2000_.massive_body(*ephemeris_, "Moon"))),
        lunar_frame_(ephemeris_.get(), moon_),
        selenocentre_(Selenocentric::origin, Velocity<Selenocentric>{}) {
    google::LogToStderr();
  }

  enum class LunarTag { surface, selenocentric };

  // This Moon-centred, Moon-fixed reference frame has the x axis pointing
  // towards the Earth, and the y axis in the direction of the velocity of the
  // Earth, see figure 1. of Russell and Lara (2006).
  using LunarSurface = Frame<LunarTag,
                             LunarTag::surface,
                             /*frame_is_inertial=*/false>;

  // This reference frame is non-rotating, with its origin at the selenocentre.
  // The axes are those of LunarSurface at J2000.
  // Note that this frame is not actually inertial, but we want to use it with
  // |KeplerOrbit|.  Perhaps we should have a concept of non-rotating, and
  // |KeplerOrbit| should check that; this is good enough for a test.
  using Selenocentric = Frame<LunarTag,
                              LunarTag::selenocentric,
                              /*frame_is_inertial=*/true>;

  // We do not use a |BodyCentredNonRotatingDynamicFrame| since that would use
  // ICRS axes.
  RigidMotion<ICRS, Selenocentric> ToSelenocentric(Instant const& t) {
    return RigidMotion<ICRS, Selenocentric>(
        RigidTransformation<ICRS, Selenocentric>(
            ephemeris_->trajectory(moon_)->EvaluatePosition(t),
            Selenocentric::origin,
            OrthogonalMap<LunarSurface, Selenocentric>::Identity() *
                lunar_frame_.ToThisFrameAtTime(J2000).orthogonal_map()),
        /*angular_velocity_of_to_frame=*/AngularVelocity<ICRS>{},
        /*velocity_of_to_frame_origin=*/
        ephemeris_->trajectory(moon_)->EvaluateVelocity(t));
  }

  SolarSystem<ICRS> const solar_system_2000_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<OblateBody<ICRS> const*> const moon_;

  BodySurfaceDynamicFrame<ICRS, LunarSurface> const lunar_frame_;
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
         /*max_degree=*/50,
         /*zonal_only=*/false,
         /*first_period_eccentricity_vector_drift=*/0.00018,
         /*first_period_descending_nodes=*/{-0.0055, +0.0051, -0.027, -0.018},
         /*period_ends=*/{+0.0026, +0.0037, -0.021, -0.0200},
         /*periods=*/10,
     },
     {
#endif
         /*max_degree=*/30,
         /*zonal_only=*/false,
         /*first_period_eccentricity_vector_drift=*/0.00032,
         /*first_period_descending_nodes=*/{-0.0058, +0.0048, -0.027, -0.018},
         /*period_ends=*/{+0.0019, +0.0050, -0.022, -0.0190},
         /*periods=*/28,
     },
     {
         /*max_degree=*/25,
         /*zonal_only=*/false,
         /*first_period_eccentricity_vector_drift=*/0.00110,
         /*first_period_descending_nodes=*/{-0.0060, +0.0044, -0.027, -0.018},
         /*period_ends=*/{-0.0017, +0.0089, -0.021, -0.0110},
         /*periods=*/28,
     },
     {
         /*max_degree=*/20,
         /*zonal_only=*/false,
         /*first_period_eccentricity_vector_drift=*/0.00130,
         /*first_period_descending_nodes=*/{-0.0064, +0.0045, -0.028, -0.018},
         /*period_ends=*/{-0.0030, +0.0100, -0.021, -0.0083},
         /*periods=*/28,
     },
     {
         /*max_degree=*/10,
         /*zonal_only=*/false,
         /*first_period_eccentricity_vector_drift=*/0.00370,
         /*first_period_descending_nodes=*/{-0.0091, +0.0036, -0.028, -0.018},
         /*period_ends=*/{-0.0160, +0.0210, -0.021, +0.0160},
         /*periods=*/28,
#if PRINCIPIA_GEOPOTENTIAL_MAX_DEGREE_50
     },
     {
         /*max_degree=*/50,
         /*zonal_only=*/true,
         /*first_period_eccentricity_vector_drift=*/0.00098,
         /*first_period_descending_nodes=*/{+0.0038, +0.0040, -0.022, -0.021},
         /*period_ends=*/{-0.0047, +0.0040, -0.025, -0.0170},
         /*periods=*/28,
#endif
     }},
};

INSTANTIATE_TEST_CASE_P(
    TruncatedSelenopotentials,
    LunarOrbitTest,
    ::testing::ValuesIn(geopotential_truncations));

TEST_P(LunarOrbitTest, NearCircularRepeatGroundTrackOrbit) {
  Time const integration_step = 10 * Second;
  LOG(INFO) << "Using a " << GetParam() << " selenopotential field";

  OFStream file(SOLUTION_DIR / "mathematica" /
                absl::StrCat("lunar_orbit_",
                             GetParam().DegreeAndOrder(),
                             ".generated.wl"));

  // We work with orbit C from Russell and Lara (2006), Repeat Ground Track
  // Lunar Orbits in the Full-Potential Plus Third-Body Problem.

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

  // The _rl values are the ones from table 1 of Russell and Lara (2006).
  Length const LU_rl = 384'400 * Kilo(Metre);
  Time const TU_rl = 375'190.258663027 * Second;
  GravitationalParameter const GM_rl =
      4'902.801076 * (Pow<3>(Kilo(Metre)) / Pow<2>(Second));

  Time const TU = 1 * Radian / moon_->angular_velocity().Norm();
  Length const LU = Cbrt((moon_->gravitational_parameter() * Pow<2>(TU)) /
                         (GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl))));
  EXPECT_THAT(moon_->gravitational_parameter() / (Pow<3>(LU) / Pow<2>(TU)),
              AlmostEquals(GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl)), 1));
  EXPECT_THAT(RelativeError(TU, TU_rl), IsNear(1.4e-3_⑴));
  EXPECT_THAT(RelativeError(LU, LU_rl), IsNear(9.0e-4_⑴));

  file << mathematica::Assign("tu", TU / Second);
  file << mathematica::Assign("lu", LU / Metre);

  Time const period = 2 * π * TU;
  int const orbits_per_period = 328;

  // Initial conditions and elements from table 2 of Russell and Lara (2006).
  Length const x0 = -4.498948742093e-03 * LU;
  Length const y0 = -1.731769313131e-03 * LU;
  Length const z0 =  0 * LU;
  Speed const  u0 = -6.203996010078e-02 * (LU / TU);
  Speed const  v0 =  7.000280770869e-02 * (LU / TU);
  Speed const  w0 =  1.588813067177e+00 * (LU / TU);

  DegreesOfFreedom<LunarSurface> const lunar_initial_state = {
      LunarSurface::origin + Displacement<LunarSurface>({x0, y0, z0}),
      Velocity<LunarSurface>({u0, v0, w0})};

  ephemeris_->Prolong(J2000);
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
                IsNear(9.0e-4_⑴));
    EXPECT_THAT(RelativeError(*initial_osculating.eccentricity, e0),
                IsNear(1.4e-10_⑴));
    EXPECT_THAT(RelativeError(initial_osculating.inclination, i0),
                IsNear(9.7e-9_⑴));
    EXPECT_THAT(RelativeError(*initial_osculating.argument_of_periapsis,
                              2 * π * Radian + ω0),
                IsNear(2.0e-11_⑴));
    EXPECT_THAT(RelativeError(initial_osculating.longitude_of_ascending_node,
                              2 * π * Radian + Ω0),
                IsNear(4.7e-13_⑴));
  }

  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(J2000, initial_state);
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRS>>(),
          integration_step));

  ephemeris_->FlowWithFixedStep(J2000 + GetParam().periods * period, *instance);

  // To find the nodes, we need to convert the trajectory to a reference frame
  // whose xy plane is the Moon's equator.
  DiscreteTrajectory<LunarSurface> surface_trajectory;
  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    surface_trajectory.Append(
        it.time(),
        lunar_frame_.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
  }

  // Drop the units when logging to Mathematica, because it is ridiculously
  // slow at parsing them.
  std::vector<double> mma_times;
  std::vector<double> mma_semimajor_axes;
  std::vector<double> mma_eccentricities;
  std::vector<double> mma_inclinations;
  std::vector<double> mma_arguments_of_periapsides;
  std::vector<double> mma_longitudes_of_ascending_nodes;
  std::vector<Vector<double, LunarSurface>> mma_displacements;

  for (Instant t = J2000; t <= J2000 + 2 * period; t += period / 50'000) {
    auto const elements = KeplerOrbit<Selenocentric>(
        *moon_,
        satellite_,
        ToSelenocentric(t)(
            trajectory.EvaluateDegreesOfFreedom(t)) - selenocentre_,
        t).elements_at_epoch();

    mma_times.push_back((t - J2000) / Second);
    mma_semimajor_axes.push_back(*elements.semimajor_axis / Metre);
    mma_inclinations.push_back(elements.inclination / Radian);
    mma_eccentricities.push_back(*elements.eccentricity);
    mma_arguments_of_periapsides.push_back(*elements.argument_of_periapsis /
                                           Radian);
    mma_longitudes_of_ascending_nodes.push_back(
        elements.longitude_of_ascending_node / Radian);
    mma_displacements.push_back(
        (surface_trajectory.EvaluatePosition(t) - LunarSurface::origin) /
        Metre);
  }

  file << mathematica::Assign("times", mma_times);
  file << mathematica::Assign("semimajorAxes", mma_semimajor_axes);
  file << mathematica::Assign("eccentricities", mma_eccentricities);
  file << mathematica::Assign("inclinations", mma_inclinations);
  file << mathematica::Assign("arguments", mma_arguments_of_periapsides);
  file << mathematica::Assign("longitudesOfAscendingNodes",
                              mma_longitudes_of_ascending_nodes);
  file << mathematica::Assign("displacements", mma_displacements);

  DiscreteTrajectory<LunarSurface> ascending_nodes;
  DiscreteTrajectory<LunarSurface> descending_nodes;
  ComputeNodes(surface_trajectory.Begin(),
               surface_trajectory.End(),
               /*north=*/Vector<double, LunarSurface>({0, 0, 1}),
               /*max_points=*/std::numeric_limits<int>::max(),
               ascending_nodes,
               descending_nodes);

  DiscreteTrajectory<ICRS> apoapsides;
  DiscreteTrajectory<ICRS> periapsides;
  ComputeApsides(*ephemeris_->trajectory(moon_),
                 trajectory.Begin(),
                 trajectory.End(),
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
    std::vector<double> mma_node_times;
    std::vector<Vector<double, LunarSurface>> mma_node_displacements;
    std::vector<double> mma_node_arguments_of_periapsides;
    std::vector<double> mma_node_eccentricities;

    for (auto it = nodes.trajectory.Begin();
         it != nodes.trajectory.End();
         ++it) {
      auto const t = it.time();
      auto const elements = KeplerOrbit<Selenocentric>(
          *moon_,
          satellite_,
          ToSelenocentric(t)(
              trajectory.EvaluateDegreesOfFreedom(t)) - selenocentre_,
          t).elements_at_epoch();

      mma_node_times.push_back((t - J2000) / Second);
      mma_node_displacements.push_back(
          (it.degrees_of_freedom().position() - LunarSurface::origin) / Metre);
      mma_node_arguments_of_periapsides.push_back(
          *elements.argument_of_periapsis / Radian);
      mma_node_eccentricities.push_back(*elements.eccentricity);

      if (nodes.name == "descending") {
        descending_node_eccentricities.push_back(*elements.eccentricity);
        descending_node_arguments.push_back(*elements.argument_of_periapsis);
      }
    }
    file << mathematica::Assign(absl::StrCat(nodes.name, "NodeTimes"),
                                mma_node_times);
    file << mathematica::Assign(absl::StrCat(nodes.name, "NodeDisplacements"),
                                mma_node_displacements);
    file << mathematica::Assign(absl::StrCat(nodes.name, "NodeArguments"),
                                mma_node_arguments_of_periapsides);
    file << mathematica::Assign(absl::StrCat(nodes.name, "NodeEccentricities"),
                                mma_node_eccentricities);
  }

  for (auto const& apsides : {Apsides{"apoapsis", apoapsides},
                              Apsides{"periapsis", periapsides}}) {
    std::vector<double> mma_apsis_times;
    std::vector<Vector<double, LunarSurface>> mma_apsis_displacements;

    for (auto it = apsides.trajectory.Begin();
         it != apsides.trajectory.End();
         ++it) {
      auto const t = it.time();

      mma_apsis_times.push_back((t - J2000) / Second);
      mma_apsis_displacements.push_back(
          (lunar_frame_.ToThisFrameAtTime(t).rigid_transformation()(
               it.degrees_of_freedom().position()) -
           LunarSurface::origin) / Metre);
    }
    file << mathematica::Assign(absl::StrCat(apsides.name, "Times"),
                                mma_apsis_times);
    file << mathematica::Assign(absl::StrCat(apsides.name, "Displacements"),
                                mma_apsis_displacements);
  }

  {
    auto const e0 = descending_node_eccentricities[0];
    auto const e1 = descending_node_eccentricities[orbits_per_period];
    auto const ω0 = descending_node_arguments[0];
    auto const ω1 = descending_node_arguments[orbits_per_period];
    EXPECT_THAT(RelativeError(GetParam().first_period_eccentricity_vector_drift,
                              Sqrt(Pow<2>(e1 * Cos(ω1) - e0 * Cos(ω0)) +
                                   Pow<2>(e1 * Sin(ω1) - e0 * Sin(ω0)))),
                IsNear(0.01_⑴));
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
        IsNear(0.01_⑴));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.max_e_cos_ω,
                      actual_first_period_descending_nodes.max_e_cos_ω),
        IsNear(0.01_⑴));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.min_e_sin_ω,
                      actual_first_period_descending_nodes.min_e_sin_ω),
        IsNear(0.01_⑴));
    EXPECT_THAT(
        RelativeError(GetParam().first_period_descending_nodes.max_e_sin_ω,
                      actual_first_period_descending_nodes.max_e_sin_ω),
        IsNear(0.01_⑴));
  }

  {
    EccentricityVectorRange actual_period_ends;
    for (int orbit = 0;
         orbit < descending_nodes.Size();
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
                IsNear(0.01_⑴));
    EXPECT_THAT(RelativeError(GetParam().period_ends.max_e_cos_ω,
                              actual_period_ends.max_e_cos_ω),
                IsNear(0.01_⑴));
    EXPECT_THAT(RelativeError(GetParam().period_ends.min_e_sin_ω,
                              actual_period_ends.min_e_sin_ω),
                IsNear(0.01_⑴));
    EXPECT_THAT(RelativeError(GetParam().period_ends.max_e_sin_ω,
                              actual_period_ends.max_e_sin_ω),
                IsNear(0.01_⑴));
  }
}

#endif

}  // namespace astronomy
}  // namespace principia
