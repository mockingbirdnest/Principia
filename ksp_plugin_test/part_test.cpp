#include "ksp_plugin/part.hpp"

#include "astronomy/epoch.hpp"
#include "base/algebra.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/instant.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/euler_solver.hpp"
#include "physics/mechanical_system.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/tensors.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace ksp_plugin {

using ::testing::_;
using ::testing::MockFunction;
using namespace principia::astronomy::_epoch;
using namespace principia::base::_algebra;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_permutation;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_identification;
using namespace principia::ksp_plugin::_part;
using namespace principia::ksp_plugin::_pile_up;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_euler_solver;
using namespace principia::physics::_mechanical_system;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_tensors;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

class PartTest : public testing::Test {
 protected:
  PartTest()
      : inertia_tensor_(MakeWaterSphereInertiaTensor(mass_)),
        part_(part_id_,
              "part",
              mass_,
              EccentricPart::origin,
              inertia_tensor_,
              RigidMotion<EccentricPart, Barycentric>::MakeNonRotatingMotion(
                  degrees_of_freedom_),
              /*deletion_callback=*/nullptr) {
    part_.apply_intrinsic_force(intrinsic_force_);
    part_.AppendToHistory(
        J2000,
        {Barycentric::origin +
             Displacement<Barycentric>({11 * Metre, 22 * Metre, 33 * Metre}),
         Velocity<Barycentric>(
             {44 * Metre / Second, 55 * Metre / Second, 66 * Metre / Second})});
  }

  DegreesOfFreedom<Barycentric> const degrees_of_freedom_ = {
      Barycentric::origin +
          Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
      Velocity<Barycentric>(
          {4 * Metre / Second, 5 * Metre / Second, 6 * Metre / Second})};
  PartId const part_id_ = 666;
  Mass const mass_ = 7 * Kilogram;
  Vector<Force, Barycentric> const intrinsic_force_ =
      Vector<Force, Barycentric>({8 * Newton, 9 * Newton, 10 * Newton});
  InertiaTensor<RigidPart> inertia_tensor_;
  Part part_;
};

// Verify that we can produce a torque that approximates the exponential decay
// of the angular velocity that PhysX implements based on `angular_drag`.  This
// behaviour makes no physical sense, but it is our duty to replicate it.
TEST(PartTestWithoutFixture, PhysXExponentialDecay) {
  // We need a non-rotating frame to construct a solver, but the frames that
  // come from the game like `RigidPart` are left-handed.  Let's have a
  // left-handed inertial frame just for this test.
  using AliceBarycentric = Frame<serialization::Frame::TestTag,
                                 Inertial,
                                 Handedness::Left,
                                 serialization::Frame::TEST>;

  // The numbers below are arbitrary, but their order of magnitude corresponds
  // to what was observed in game, so as to give us an idea of the "quality" of
  // our approximation.
  Instant const t0;
  Time const Δt = 20 * Milli(Second);

  InertiaTensor<RigidPart> const inertia_tensor(
      {{
           10.0e3 * si::Unit<MomentOfInertia>,
           -5.5e3 * si::Unit<MomentOfInertia>,
           2.0e3 * si::Unit<MomentOfInertia>,
       },
       {
           -5.5e3 * si::Unit<MomentOfInertia>,
           14.0e3 * si::Unit<MomentOfInertia>,
           -6.0e3 * si::Unit<MomentOfInertia>,
       },
       {
           2.0e3 * si::Unit<MomentOfInertia>,
           -6.0e3 * si::Unit<MomentOfInertia>,
           7.0e3 * si::Unit<MomentOfInertia>,
       }});

  auto const eigensystem = inertia_tensor.Diagonalize<PartPrincipalAxes>();
  EXPECT_THAT(
      eigensystem.form.coordinates().Diagonal(),
      AlmostEquals(R3Element<MomentOfInertia>(
                       {3.32030832353663583474e3 * si::Unit<MomentOfInertia>,
                        7.0799313792826743057e3 * si::Unit<MomentOfInertia>,
                        20.5997602971806898595e3 * si::Unit<MomentOfInertia>}),
                   2));

  Rotation<PartPrincipalAxes, RigidPart> const initial_attitude =
      eigensystem.rotation;

  AngularVelocity<RigidPart> const angular_velocity(
      {0.1 * Radian / Second, -0.02 * Radian / Second, 0.05 * Radian / Second});

  Inverse<Time> const angular_drag = 20 / Second;

  // This is the code being tested.
  Bivector<Torque, RigidPart> const intrinsic_torque =
      Part::DragTorqueFromAngularVelocity(
          angular_drag, Δt, angular_velocity, inertia_tensor);

  Bivector<AngularMomentum, RigidPart> const initial_angular_momentum =
      inertia_tensor * angular_velocity;
  Bivector<AngularMomentum, RigidPart> const updated_angular_momentum =
      initial_angular_momentum + intrinsic_torque * Δt;

  Rotation<RigidPart, AliceBarycentric> const rigid_part_to_alice_barycentric =
      Rotation<RigidPart, AliceBarycentric>::Identity();
  Rotation<AliceBarycentric, RigidPart> const alice_barycentric_to_rigid_part =
      rigid_part_to_alice_barycentric.Inverse();

  EulerSolver<AliceBarycentric, PartPrincipalAxes> const euler_solver(
      eigensystem.form,
      rigid_part_to_alice_barycentric(updated_angular_momentum),
      rigid_part_to_alice_barycentric * initial_attitude,
      t0);

  Bivector<AngularMomentum, PartPrincipalAxes> const angular_momentum_before =
      euler_solver.AngularMomentumAt(t0);
  Rotation<PartPrincipalAxes, RigidPart> const attitude_before =
      alice_barycentric_to_rigid_part *
      euler_solver.AttitudeAt(angular_momentum_before, t0);
  AngularVelocity<PartPrincipalAxes> const angular_velocity_before =
      euler_solver.AngularVelocityFor(angular_momentum_before);

  AngularVelocity<RigidPart> const expected_angular_velocity_before =
      angular_velocity * (1 - angular_drag * Δt);
  AngularVelocity<RigidPart> const actual_angular_velocity_before =
      attitude_before(angular_velocity_before);
  EXPECT_THAT(actual_angular_velocity_before.Norm(),
              RelativeErrorFrom(expected_angular_velocity_before.Norm(),
                                IsNear(1.1e-3_(1))));
  EXPECT_THAT(AngleBetween(actual_angular_velocity_before,
                           expected_angular_velocity_before),
              IsNear(6.4e-3_(1) * Radian));

  Bivector<AngularMomentum, PartPrincipalAxes> const angular_momentum_after =
      euler_solver.AngularMomentumAt(t0 + Δt);
  Rotation<PartPrincipalAxes, RigidPart> const attitude_after =
      alice_barycentric_to_rigid_part *
      euler_solver.AttitudeAt(angular_momentum_after, t0 + Δt);
  AngularVelocity<PartPrincipalAxes> const angular_velocity_after =
      euler_solver.AngularVelocityFor(angular_momentum_after);

  AngularVelocity<RigidPart> const expected_angular_velocity_after =
      angular_velocity * (1 - angular_drag * Δt);
  AngularVelocity<RigidPart> const actual_angular_velocity_after =
      attitude_after(angular_velocity_after);
  EXPECT_THAT(actual_angular_velocity_after.Norm(),
              RelativeErrorFrom(expected_angular_velocity_after.Norm(),
                                IsNear(7.1e-4_(1))));
  EXPECT_THAT(AngleBetween(actual_angular_velocity_after,
                           expected_angular_velocity_after),
              IsNear(4.1e-3_(1) * Radian));
}

TEST_F(PartTest, Serialization) {
  MockFunction<int(not_null<not_null<PileUp const*>>)>
      serialization_index_for_pile_up;
  EXPECT_CALL(serialization_index_for_pile_up, Call(_)).Times(0);

  serialization::Part message;
  part_.WriteToMessage(&message,
                       serialization_index_for_pile_up.AsStdFunction());
  EXPECT_EQ(part_id_, message.part_id());
  EXPECT_TRUE(message.has_inertia_tensor());
  EXPECT_TRUE(message.has_mass());
  EXPECT_EQ(7, message.mass().magnitude());
  EXPECT_TRUE(message.has_intrinsic_force());
  EXPECT_TRUE(message.intrinsic_force().has_vector());
  EXPECT_EQ(8, message.intrinsic_force().vector().x().quantity().magnitude());
  EXPECT_EQ(9, message.intrinsic_force().vector().y().quantity().magnitude());
  EXPECT_EQ(10, message.intrinsic_force().vector().z().quantity().magnitude());
  EXPECT_TRUE(message.has_rigid_motion());
  EXPECT_TRUE(message.rigid_motion().has_rigid_transformation());
  EXPECT_TRUE(message.rigid_motion().rigid_transformation().has_to_origin());
  EXPECT_TRUE(message.rigid_motion()
                  .rigid_transformation()
                  .to_origin()
                  .has_multivector());
  EXPECT_TRUE(message.rigid_motion()
                  .rigid_transformation()
                  .to_origin()
                  .multivector()
                  .has_vector());
  EXPECT_EQ(1,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .x()
                .quantity()
                .magnitude());
  EXPECT_EQ(2,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .y()
                .quantity()
                .magnitude());
  EXPECT_EQ(3,
            message.rigid_motion()
                .rigid_transformation()
                .to_origin()
                .multivector()
                .vector()
                .z()
                .quantity()
                .magnitude());
  EXPECT_TRUE(
      message.rigid_motion().velocity_of_to_frame_origin().has_vector());
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .x()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-4, 6));
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .y()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-6, 2));
  EXPECT_THAT(message.rigid_motion()
                  .velocity_of_to_frame_origin()
                  .vector()
                  .z()
                  .quantity()
                  .magnitude(),
              AlmostEquals(-5, 2));
  EXPECT_EQ(1, message.prehistory().segment_size());
  EXPECT_EQ(1, message.prehistory().segment(0).zfp().timeline_size());

  auto const p = Part::ReadFromMessage(message, /*deletion_callback=*/nullptr);
  EXPECT_EQ(part_.inertia_tensor(), p->inertia_tensor());
  EXPECT_EQ(part_.intrinsic_force(), p->intrinsic_force());
  EXPECT_EQ(part_.rigid_motion()({RigidPart::origin, RigidPart::unmoving}),
            p->rigid_motion()({RigidPart::origin, RigidPart::unmoving}));
  EXPECT_EQ(part_.rigid_motion().angular_velocity_of<RigidPart>(),
            p->rigid_motion().angular_velocity_of<RigidPart>());

  serialization::Part second_message;
  p->WriteToMessage(&second_message,
                    serialization_index_for_pile_up.AsStdFunction());
  EXPECT_THAT(message, EqualsProto(second_message));
}

}  // namespace ksp_plugin
}  // namespace principia
