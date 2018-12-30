
#include <filesystem>
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
using physics::ComputeNodes;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
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
using quantities::GravitationalParameter;
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
using testing_utilities::AlmostEquals;
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::Slope;

namespace astronomy {

class LunarOrbitTest : public ::testing::Test {
 protected:
  LunarOrbitTest()
      : solar_system_2000_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt"),
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
        instantaneous_moon_(InstantaneousLunarSurface::origin,
                            Velocity<InstantaneousLunarSurface>{}) {
    google::LogToStderr();
  }

  enum class LunarSurfaceTag { rotating, instantaneous_inertial };

  // This Moon-centred, Moon-fixed reference frame has the x axis pointing
  // towards the Earth, and the y axis in the direction of the velocity of the
  // Earth, see figure 1. of Russell and Lara (2006).
  using LunarSurface = Frame<LunarSurfaceTag,
                             LunarSurfaceTag::rotating,
                             /*frame_is_inertial=*/false>;

  // At any time t, |ToInstantaneousLunarSurfaceFrame(t)| converts to an
  // inertial frame wherein the Moon is immobile at the origin at t, the x axis
  // points towards the position of the Earth at t, and the y axis is in the
  // direction of the velocity of the Earth at t.  In other words, this frame is
  // the inertial continuation of |LunarSurface| at t.  Note that this type
  // represents a different reference frame at each time, instead of a single
  // reference frame.  Time evolutions should not be computed in these frames.
  using InstantaneousLunarSurface =
      Frame<LunarSurfaceTag,
            LunarSurfaceTag::instantaneous_inertial,
            /*frame_is_inertial=*/true>;

  RigidMotion<ICRS, InstantaneousLunarSurface>
  ToInstantaneousLunarSurfaceFrame(Instant const& t) {
    RigidTransformation<LunarSurface, InstantaneousLunarSurface>
        trivial_rigid_transformation_at_t(
            LunarSurface::origin,
            InstantaneousLunarSurface::origin,
            OrthogonalMap<LunarSurface, InstantaneousLunarSurface>::Identity());
    return RigidMotion<ICRS, InstantaneousLunarSurface>(
        trivial_rigid_transformation_at_t *
            lunar_frame_.ToThisFrameAtTime(t).rigid_transformation(),
        /*angular_velocity_of_to_frame=*/AngularVelocity<ICRS>{},
        /*velocity_of_to_frame_origin=*/
        ephemeris_->trajectory(moon_)->EvaluateVelocity(t));
  }

  SolarSystem<ICRS> const solar_system_2000_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<OblateBody<ICRS> const*> const moon_;

  BodySurfaceDynamicFrame<ICRS, LunarSurface> const lunar_frame_;
  DegreesOfFreedom<InstantaneousLunarSurface> instantaneous_moon_;

  MasslessBody const satellite_;
};

#if !defined(_DEBUG)

TEST_F(LunarOrbitTest, NearCircularRepeatGroundTrackOrbit) {
  Time const integration_step = 10 * Second;

  // We work with orbit C from Russell and Lara (2006), Repeat Ground Track Lunar
  // Orbits in the Full-Potential Plus Third-Body Problem.

  // The length and time units LU and TU are such that, in an idealized
  // Earth-Moon system, the Earth-Moon distance is 1 LU and the angular
  // frequency of the body-fixed moon frame is θ′ = 1 rad / TU, see figure 1 and
  // table 1 of Russell and Lara (2006).
  // In order to best reproduce the results of the paper, we choose our TU such
  // that the rotational period of the moon is TU, and our LU such that the
  // Moon's gravitational parameter has the same value in LU³/TU² as the
  // gravitational parameter used in the paper.
  // With cartesian initial conditions in the surface frame, these two
  // properties ensure that the initial osculating lunar orbit has the same
  // orientation, eccentricity, and anomaly.
  Length const LU_rl = 384'400 * Kilo(Metre);
  Time const TU_rl = 375'190.258663027 * Second;
  GravitationalParameter const GM_rl =
      4'902.801076 * (Pow<3>(Kilo(Metre)) / Pow<2>(Second));

  Time const TU = 1 * Radian / moon_->angular_velocity().Norm();
  Length const LU = Cbrt((moon_->gravitational_parameter() * Pow<2>(TU)) /
                         (GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl))));
  EXPECT_THAT(moon_->gravitational_parameter() / (Pow<3>(LU) / Pow<2>(TU)),
              AlmostEquals(GM_rl / (Pow<3>(LU_rl) / Pow<2>(TU_rl)), 0));
  EXPECT_THAT(RelativeError(TU, TU_rl), IsNear(1.4e-3));
  EXPECT_THAT(RelativeError(LU, LU_rl), IsNear(0));

  // Initial conditions and elements from table 2 of Russell and Lara (2006).
  Length const x0 = -4.498948742093e-03 * LU;
  Length const y0 = -1.731769313131e-03 * LU;
  Length const z0 =  0 * LU;
  Speed const u0  = -6.203996010078e-02 * (LU / TU);
  Speed const v0  =  7.000280770869e-02 * (LU / TU);
  Speed const w0  =  1.588813067177e+00 * (LU / TU);

  DegreesOfFreedom<LunarSurface> lunar_initial_state = {
      LunarSurface::origin + Displacement<LunarSurface>({x0, y0, z0}),
      Velocity<LunarSurface>({u0, v0, w0})};

  ephemeris_->Prolong(J2000);
  DegreesOfFreedom<ICRS> initial_state =
      lunar_frame_.FromThisFrameAtTime(J2000)(lunar_initial_state);

  {
    KeplerOrbit<InstantaneousLunarSurface> initial_orbit(
        *moon_,
        satellite_,
        ToInstantaneousLunarSurfaceFrame(J2000)(initial_state) -
            instantaneous_moon_,
        J2000);
    EXPECT_THAT(RelativeError(*initial_orbit.elements_at_epoch().semimajor_axis,
                              +1.861791339407e+03 * Kilo(Metre)),
                IsNear(2.4e-3));
    EXPECT_THAT(RelativeError(*initial_orbit.elements_at_epoch().eccentricity,
                              +2.110475283361e-02),
                IsNear(1.9e-2));
    EXPECT_THAT(AbsoluteError(initial_orbit.elements_at_epoch().inclination,
                              +9.298309294740e+01 * Degree),
                IsNear(10 * ArcMinute));
    EXPECT_THAT(
        AbsoluteError(*initial_orbit.elements_at_epoch().argument_of_periapsis,
                      7.839337618501e+01 * Degree),
        IsNear(6.5 * Degree));
    EXPECT_THAT(
        RelativeError(
            initial_orbit.elements_at_epoch().longitude_of_ascending_node,
            2 * π * Radian - 1.589469097527e+02 * Degree),
        IsNear(2.4e-4));
  }

  Time const integration_duration = 2 * 28 * Day;

  DiscreteTrajectory<ICRS> trajectory;
  trajectory.Append(J2000, initial_state);
  auto const instance = ephemeris_->NewInstance(
      {&trajectory},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRS>>(),
          integration_step));

  // Remember that because of #228 we need to loop over FlowWithFixedStep.
  ephemeris_->FlowWithFixedStep(J2000 + integration_duration, *instance);

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
  base::OFStream file(SOLUTION_DIR / "mathematica" /
                      "lunar_orbit.generated.wl");
  std::vector<Vector<double, LunarSurface>> mma_displacements;
  std::vector<double> mma_arguments_of_periapsides;
  std::vector<double> mma_eccentricities;
  std::vector<double> mma_times;

  for (Instant t = J2000; t <= J2000 + integration_duration;
       t += integration_duration / 100'000.0) {
    auto const elements = KeplerOrbit<InstantaneousLunarSurface>(
        *moon_,
        satellite_,
        ToInstantaneousLunarSurfaceFrame(t)(
            trajectory.EvaluateDegreesOfFreedom(t)) - instantaneous_moon_,
        t).elements_at_epoch();

    mma_times.push_back((t - J2000) / Second);
    mma_displacements.push_back(
        (surface_trajectory.EvaluatePosition(t) - LunarSurface::origin) /
        Metre);
    mma_arguments_of_periapsides.push_back(*elements.argument_of_periapsis /
                                           Radian);
    mma_eccentricities.push_back(*elements.eccentricity);
  }

  file << mathematica::Assign("times", mma_times);
  file << mathematica::Assign("displacements", mma_displacements);
  file << mathematica::Assign("arguments", mma_arguments_of_periapsides);
  file << mathematica::Assign("eccentricities", mma_eccentricities);

  DiscreteTrajectory<LunarSurface> ascending_nodes;
  DiscreteTrajectory<LunarSurface> descending_nodes;
  ComputeNodes(surface_trajectory.Begin(),
               surface_trajectory.End(),
               /*north=*/Vector<double, LunarSurface>({0, 0, 1}),
               ascending_nodes,
               descending_nodes);

  struct Nodes {
    std::string_view const name;
    DiscreteTrajectory<LunarSurface> const& trajectory;
  };

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
      auto const elements = KeplerOrbit<InstantaneousLunarSurface>(
          *moon_,
          satellite_,
          ToInstantaneousLunarSurfaceFrame(t)(
              trajectory.EvaluateDegreesOfFreedom(t)) - instantaneous_moon_,
          t).elements_at_epoch();

      mma_node_times.push_back((t - J2000) / Second);
      mma_node_displacements.push_back(
          (it.degrees_of_freedom().position() - LunarSurface::origin) / Metre);
      mma_node_arguments_of_periapsides.push_back(
          *elements.argument_of_periapsis / Radian);
      mma_node_eccentricities.push_back(*elements.eccentricity);
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
}

#endif

}  // namespace astronomy
}  // namespace principia
