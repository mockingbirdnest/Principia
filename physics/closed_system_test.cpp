
#include "physics/closed_system.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/componentwise.hpp"

namespace principia {

using geometry::Bivector;
using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
using geometry::OrthogonalMap;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::Vector;
using geometry::Velocity;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::AngularMomentum;
using quantities::Energy;
using quantities::Length;
using quantities::Mass;
using quantities::MomentOfInertia;
using quantities::Momentum;
using quantities::Pow;
using quantities::Speed;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::Componentwise;
using ::testing::Eq;

namespace physics {

class ClosedSystemTest : public testing::Test{
 protected:
  using InertialFrame = Frame<enum class InertialTag, Inertial>;
  using SystemFrame = Frame<enum class SystemTag>;

  ClosedSystem<InertialFrame, SystemFrame> system_;
};

TEST_F(ClosedSystemTest, TwoPointMasses) {
  constexpr Mass m1 = 1 * Kilogram;
  constexpr Mass m2 = 7 * Kilogram;
  constexpr Length r = 5 * Metre;
  constexpr Speed v = 3 * Metre / Second;

  using M1 = Frame<enum class M1Tag>;
  RigidMotion<M1, InertialFrame> m1_motion(
      RigidTransformation<M1, InertialFrame>(
          M1::origin,
          InertialFrame::origin,
          OrthogonalMap<M1, InertialFrame>::Identity()),
      InertialFrame::nonrotating,
      InertialFrame::unmoving);
  system_.AddRigidBody(
      m1_motion, m1, SymmetricBilinearForm<MomentOfInertia, M1>{});

  using M2 = Frame<enum class M2Tag>;
  RigidMotion<M2, InertialFrame> m2_motion(
      RigidTransformation<M2, InertialFrame>(
          M2::origin,
          InertialFrame::origin +
              Displacement<InertialFrame>({r, 0 * Metre, 0 * Metre}),
          OrthogonalMap<M2, InertialFrame>::Identity()),
      InertialFrame::nonrotating,
      Velocity<InertialFrame>({0 * Metre / Second, v, 0 * Metre / Second}));
  system_.AddRigidBody(
      m2_motion, m2, SymmetricBilinearForm<MomentOfInertia, M2>{});

  constexpr Mass μ = m1 * m2 / (m1 + m2);

  EXPECT_THAT(system_.mass(), Eq(m1 + m2));
  // TODO(egg): this would be useful in general.
  DegreesOfFreedom<InertialFrame> const unmoving_origin = {
      InertialFrame::origin, InertialFrame::unmoving};

  constexpr double ⅞ = 7.0 / 8;
  EXPECT_THAT(
      system_.linear_motion()({SystemFrame::origin, SystemFrame::unmoving}) -
          unmoving_origin,
      Componentwise(
          Componentwise(⅞ * 5 * Metre, 0 * Metre, 0 * Metre),
          Componentwise(
              0 * Metre / Second, ⅞ * 3 * Metre / Second, 0 * Metre / Second)));
  EXPECT_THAT(
      system_.angular_momentum(),
      Componentwise(AngularMomentum{}, AngularMomentum{}, r * μ * v * Radian));
  EXPECT_THAT(system_.inertia_tensor().coordinates(),
              Eq(R3x3Matrix<MomentOfInertia>({μ * Pow<2>(r), {}, {}}, {}, {})));

  Mass const m = system_.mass();
  Vector<Momentum, InertialFrame> const p =
      m * system_.linear_motion()({SystemFrame::origin, SystemFrame::unmoving})
              .velocity();
  Bivector<AngularMomentum, SystemFrame> const L = system_.angular_momentum();
  Bivector<double, SystemFrame> const axis = Normalize(L);
  SymmetricBilinearForm<MomentOfInertia, SystemFrame> I =
      system_.inertia_tensor();
  // Compute the kinetic energy of the system in |InertialFrame| as the energy
  // of its linear motion plus that of its rotational motion, and compare that
  // with sum of the kinetic energies of the point masses.  Note that m1 is
  // stationary in |InertialFrame|.
  EXPECT_THAT(p.Norm²() / (2 * m) +
                  (L / Radian).Norm²() /
                      (2 * InnerProduct(axis, Anticommutator(I, axis))),
              Eq(0.5 * m2 * Pow<2>(v)));
}

}  // namespace physics
}  // namespace principia
