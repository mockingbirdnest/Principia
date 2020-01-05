
#include "physics/closed_system.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "testing_utilities/componentwise.hpp"

namespace principia {

using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
using geometry::OrthogonalMap;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::Velocity;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::AngularMomentum;
using quantities::Length;
using quantities::Mass;
using quantities::MomentOfInertia;
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

  EXPECT_THAT(
      system_.angular_momentum(),
      Componentwise(AngularMomentum{}, AngularMomentum{}, r * μ * v * Radian));
  EXPECT_THAT(system_.inertia_tensor().coordinates(),
              Eq(R3x3Matrix<MomentOfInertia>({μ * Pow<2>(r), {}, {}}, {}, {})));
}

}  // namespace physics
}  // namespace principia
