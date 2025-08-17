#include "physics/mechanical_system.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace physics {

using ::testing::Eq;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_mechanical_system;
using namespace principia::physics::_rigid_motion;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_componentwise;

class MechanicalSystemTest : public testing::Test{
 protected:
  using InertialFrame = Frame<struct InertialTag, Inertial>;
  using SystemFrame = Frame<struct SystemTag, NonRotating>;

  MechanicalSystem<InertialFrame, SystemFrame> system_;
};

TEST_F(MechanicalSystemTest, RigidTwoPointMasses) {
  // A rigid system of two point masses.
  constexpr Mass m1 = 1 * Kilogram;
  constexpr Mass m2 = 7 * Kilogram;
  constexpr Length r = 5 * Metre;
  constexpr Speed v = 3 * Metre / Second;

  using M1 = Frame<struct M1Tag>;
  RigidMotion<M1, InertialFrame> m1_motion(
      RigidTransformation<M1, InertialFrame>(
          M1::origin,
          InertialFrame::origin,
          OrthogonalMap<M1, InertialFrame>::Identity()),
      InertialFrame::nonrotating,
      InertialFrame::unmoving);
  system_.AddRigidBody(
      m1_motion, m1, SymmetricBilinearForm<MomentOfInertia, M1, Bivector>{});

  using M2 = Frame<struct M2Tag>;
  RigidMotion<M2, InertialFrame> m2_motion(
      RigidTransformation<M2, InertialFrame>(
          M2::origin,
          InertialFrame::origin +
              Displacement<InertialFrame>({r, 0 * Metre, 0 * Metre}),
          OrthogonalMap<M2, InertialFrame>::Identity()),
      InertialFrame::nonrotating,
      Velocity<InertialFrame>({0 * Metre / Second, v, 0 * Metre / Second}));
  system_.AddRigidBody(
      m2_motion, m2, SymmetricBilinearForm<MomentOfInertia, M2, Bivector>{});

  constexpr Mass μ = m1 * m2 / (m1 + m2);

  EXPECT_THAT(system_.mass(), Eq(m1 + m2));
  // TODO(egg): this would be useful in general.
  DegreesOfFreedom<InertialFrame> const unmoving_origin = {
      InertialFrame::origin, InertialFrame::unmoving};

  EXPECT_THAT(
      system_.LinearMotion()({SystemFrame::origin, SystemFrame::unmoving}) -
          unmoving_origin,
      Componentwise(
          Componentwise(7.0 * r / 8.0, 0 * Metre, 0 * Metre),
          Componentwise(
              0 * Metre / Second, 7.0 * v / 8.0, 0 * Metre / Second)));
  EXPECT_THAT(
      system_.AngularMomentum(),
      Componentwise(AngularMomentum{}, AngularMomentum{}, r * μ * v * Radian));
  EXPECT_THAT(system_.InertiaTensor().coordinates(),
              Eq(R3x3Matrix<MomentOfInertia>::DiagonalMatrix(
                  {{}, μ * Pow<2>(r), μ * Pow<2>(r)})));

  Mass const m = system_.mass();
  Vector<Momentum, InertialFrame> const p =
      m * system_.LinearMotion()({SystemFrame::origin, SystemFrame::unmoving})
              .velocity();
  Bivector<AngularMomentum, SystemFrame> const L = system_.AngularMomentum();
  Bivector<double, SystemFrame> const axis = Normalize(L);
  SymmetricBilinearForm<MomentOfInertia, SystemFrame, Bivector> const I =
      system_.InertiaTensor();
  // Compute the kinetic energy of the system in `InertialFrame` as the energy
  // of its linear motion plus that of its rotational motion, and compare that
  // with sum of the kinetic energies of the point masses.  Note that m1 is
  // stationary in `InertialFrame`.
  EXPECT_THAT(p.Norm²() / (2 * m) + (L / Radian).Norm²() / (2 * I(axis, axis)),
              Eq(0.5 * m2 * Pow<2>(v)));
}

TEST_F(MechanicalSystemTest, RigidTwoCubes) {
  // A rigid system of two lead cubes with three-metre-long sides, joined in a
  // cuboid 3 m × 3 m × 6 m.
  constexpr Mass cube_mass = 10 * Tonne;
  constexpr Length cube_side = 3 * Metre;
  using Cube = Frame<struct CubeTag>;
  SymmetricBilinearForm<MomentOfInertia, Cube, Bivector> const cube_inertia(
      R3x3Matrix<MomentOfInertia>::DiagonalMatrix(
          {cube_mass * Pow<2>(cube_side) / 6,
           cube_mass * Pow<2>(cube_side) / 6,
           cube_mass * Pow<2>(cube_side) / 6}));
  AngularVelocity<InertialFrame> ω(
      {0 * Radian / Second, 0 * Radian / Second, 1 * Radian / Second});

  system_.AddRigidBody(RigidMotion<Cube, InertialFrame>(
                           RigidTransformation<Cube, InertialFrame>(
                               Cube::origin,
                               InertialFrame::origin +
                                   Displacement<InertialFrame>(
                                       {-cube_side / 2, 0 * Metre, 0 * Metre}),
                               OrthogonalMap<Cube, InertialFrame>::Identity()),
                           ω,
                           Velocity<InertialFrame>({0 * Metre / Second,
                                                    -cube_side / 2 / Second,
                                                    0 * Metre / Second})),
                       cube_mass,
                       cube_inertia);
  system_.AddRigidBody(RigidMotion<Cube, InertialFrame>(
                           RigidTransformation<Cube, InertialFrame>(
                               Cube::origin,
                               InertialFrame::origin +
                                   Displacement<InertialFrame>(
                                       {cube_side / 2, 0 * Metre, 0 * Metre}),
                               OrthogonalMap<Cube, InertialFrame>::Identity()),
                           ω,
                           Velocity<InertialFrame>({0 * Metre / Second,
                                                    cube_side / 2 / Second,
                                                    0 * Metre / Second})),
                       cube_mass,
                       cube_inertia);

  EXPECT_THAT(system_.mass(), Eq(2 * cube_mass));
  EXPECT_THAT(
      system_.InertiaTensor().coordinates(),
      Eq(R3x3Matrix<MomentOfInertia>::DiagonalMatrix(
          {2 * cube_mass * (2 * Pow<2>(cube_side)) / 12,
           2 * cube_mass * (Pow<2>(2 * cube_side) + Pow<2>(cube_side)) / 12,
           2 * cube_mass * (Pow<2>(2 * cube_side) + Pow<2>(cube_side)) / 12})));
  EXPECT_THAT(
      system_.AngularMomentum(),
      Eq(system_.InertiaTensor() * Identity<InertialFrame, SystemFrame>()(ω)));
}

TEST_F(MechanicalSystemTest, NonRigidTwoCubes) {
  // A system of two counter-rotating lead cubes.
  constexpr Mass cube_mass = 10 * Tonne;
  constexpr Length cube_side = 3 * Metre;
  using Cube = Frame<struct CubeTag>;
  SymmetricBilinearForm<MomentOfInertia, Cube, Bivector> const cube_inertia(
      R3x3Matrix<MomentOfInertia>::DiagonalMatrix(
          {cube_mass * Pow<2>(cube_side) / 6,
           cube_mass * Pow<2>(cube_side) / 6,
           cube_mass * Pow<2>(cube_side) / 6}));
  AngularVelocity<InertialFrame> ω(
      {1 * Radian / Second, 0 * Radian / Second, 0 * Radian / Second});

  system_.AddRigidBody(RigidMotion<Cube, InertialFrame>(
                           RigidTransformation<Cube, InertialFrame>(
                               Cube::origin,
                               InertialFrame::origin +
                                   Displacement<InertialFrame>(
                                       {-cube_side / 2, 0 * Metre, 0 * Metre}),
                               OrthogonalMap<Cube, InertialFrame>::Identity()),
                           ω,
                           InertialFrame::unmoving),
                       cube_mass,
                       cube_inertia);
  system_.AddRigidBody(RigidMotion<Cube, InertialFrame>(
                           RigidTransformation<Cube, InertialFrame>(
                               Cube::origin,
                               InertialFrame::origin +
                                   Displacement<InertialFrame>(
                                       {cube_side / 2, 0 * Metre, 0 * Metre}),
                               OrthogonalMap<Cube, InertialFrame>::Identity()),
                           -ω,
                           InertialFrame::unmoving),
                       cube_mass,
                       cube_inertia);

  EXPECT_THAT(system_.mass(), Eq(2 * cube_mass));
  EXPECT_THAT(
      system_.InertiaTensor().coordinates(),
      Eq(R3x3Matrix<MomentOfInertia>::DiagonalMatrix(
          {2 * cube_mass * (2 * Pow<2>(cube_side)) / 12,
           2 * cube_mass * (Pow<2>(2 * cube_side) + Pow<2>(cube_side)) / 12,
           2 * cube_mass * (Pow<2>(2 * cube_side) + Pow<2>(cube_side)) / 12})));
  EXPECT_THAT(system_.AngularMomentum(),
              Eq(Bivector<AngularMomentum, SystemFrame>{}));
}

}  // namespace physics
}  // namespace principia
