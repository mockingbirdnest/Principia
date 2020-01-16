
#include "geometry/rotation.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace geometry {
namespace internal_rotation {

using quantities::ArcCos;
using quantities::ArcSin;
using quantities::Length;
using quantities::Pow;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::Eq;
using ::testing::Gt;
using ::testing::Lt;

class RotationTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;
  using World1 = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST1>;
  using Orth = OrthogonalMap<World, World>;
  using Rot = Rotation<World, World>;

  RotationTest()
      : vector_(Vector<Length, World>(
            R3Element<Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        bivector_(Bivector<Length, World>(
            R3Element<Length>(
                1.0 * Metre, 2.0 * Metre, 3.0 * Metre))),
        trivector_(Trivector<Length, World>(4.0 * Metre)),
        form_(SymmetricBilinearForm<Length, World>(
            R3x3Matrix<Length>({1.0 * Metre, 2.0 * Metre, 3.0 * Metre},
                               {2.0 * Metre, -5.0 * Metre, 6.0 * Metre},
                               {3.0 * Metre, 6.0 * Metre, 4.0 * Metre}))),
        e1_(Vector<double, World>(R3Element<double>({1, 0, 0}))),
        e2_(Vector<double, World>(R3Element<double>({0, 1, 0}))),
        e3_(Vector<double, World>(R3Element<double>({0, 0, 1}))),
        rotation_a_(Rot(120 * Degree, Bivector<double, World>({1, 1, 1}))),
        rotation_b_(Rot(90 * Degree, Bivector<double, World>({1, 0, 0}))),
        rotation_c_(Rot(ToQuaternion(
                            R3x3Matrix<double>({{0.5, 0.5 * sqrt(3), 0},
                                                {-0.5 * sqrt(3), 0.5, 0},
                                                {0, 0, 1}})))) {}

  Vector<Length, World> const vector_;
  Bivector<Length, World> const bivector_;
  Trivector<Length, World> const trivector_;
  SymmetricBilinearForm<Length, World> const form_;
  Vector<double, World> e1_;
  Vector<double, World> e2_;
  Vector<double, World> e3_;
  Rot rotation_a_;
  Rot rotation_b_;
  Rot rotation_c_;
};

using RotationDeathTest = RotationTest;

TEST_F(RotationTest, Identity) {
  EXPECT_THAT(vector_, Eq(Rot::Identity()(vector_)));
  EXPECT_THAT(bivector_, Eq(Rot::Identity()(bivector_)));
  EXPECT_THAT(trivector_, Eq(Rot::Identity()(trivector_)));
}

TEST_F(RotationTest, AppliedToVector) {
  EXPECT_THAT(rotation_a_(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(3.0 * Metre,
                                    1.0 * Metre,
                                    2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(1.0 * Metre,
                                    -3.0 * Metre,
                                    2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>((0.5 + sqrt(3.0)) * Metre,
                                    (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                    3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToBivector) {
  EXPECT_THAT(rotation_a_(bivector_),
              AlmostEquals(Bivector<Length, World>(
                  R3Element<Length>(3.0 * Metre,
                                    1.0 * Metre,
                                    2.0 * Metre)), 4));
  EXPECT_THAT(rotation_b_(bivector_),
              AlmostEquals(Bivector<Length, World>(
                  R3Element<Length>(1.0 * Metre,
                                    -3.0 * Metre,
                                    2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_(bivector_),
              AlmostEquals(Bivector<Length, World>(
                  R3Element<Length>((0.5 + sqrt(3.0)) * Metre,
                                    (1.0 - 0.5 * sqrt(3.0)) * Metre,
                                    3.0 * Metre)), 0));
}

TEST_F(RotationTest, AppliedToTrivector) {
  EXPECT_THAT(rotation_a_(trivector_),
              AlmostEquals(Trivector<Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_b_(trivector_),
              AlmostEquals(Trivector<Length, World>(
                  4.0 * Metre), 0));
  EXPECT_THAT(rotation_c_(trivector_),
              AlmostEquals(Trivector<Length, World>(
                  4.0 * Metre), 0));
}

TEST_F(RotationTest, AppliedToSymmetricBilinearForm) {
  EXPECT_THAT(form_(vector_, vector_),
              AlmostEquals(115.0 * Pow<3>(Metre), 0));
  EXPECT_THAT(rotation_a_(form_)(rotation_a_(vector_), rotation_a_(vector_)),
              AlmostEquals(115.0 * Pow<3>(Metre), 0, 1));
  EXPECT_THAT(rotation_b_(form_)(rotation_b_(vector_), rotation_b_(vector_)),
              AlmostEquals(115.0 * Pow<3>(Metre), 1, 5));
  EXPECT_THAT(rotation_c_(form_)(rotation_c_(vector_), rotation_c_(vector_)),
              AlmostEquals(115.0 * Pow<3>(Metre), 0));
}

TEST_F(RotationTest, Determinant) {
  EXPECT_TRUE(rotation_a_.Determinant().is_positive());
  EXPECT_TRUE(rotation_b_.Determinant().is_positive());
  EXPECT_TRUE(rotation_c_.Determinant().is_positive());
}

TEST_F(RotationTest, Inverse) {
  EXPECT_THAT(rotation_a_.Inverse()(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(2.0 * Metre,
                                    3.0 * Metre,
                                    1.0 * Metre)), 2));
  EXPECT_THAT(rotation_b_.Inverse()(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(1.0 * Metre,
                                    3.0 * Metre,
                                    -2.0 * Metre)), 1, 2));
  EXPECT_THAT(rotation_c_.Inverse()(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>((0.5 - sqrt(3.0)) * Metre,
                                    (1.0 + 0.5 * sqrt(3.0)) * Metre,
                                    3.0 * Metre)), 0));
}

TEST_F(RotationTest, Composition) {
  Rot const rotation_ab = rotation_a_ * rotation_b_;
  EXPECT_THAT(rotation_ab(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(2.0 * Metre,
                                                1.0 * Metre,
                                                -3.0 * Metre)), 4, 6));
}

TEST_F(RotationTest, Forget) {
  Orth const orthogonal_a = rotation_a_.Forget<OrthogonalMap>();
  EXPECT_THAT(orthogonal_a(vector_),
              AlmostEquals(Vector<Length, World>(
                  R3Element<Length>(3.0 * Metre,
                                                1.0 * Metre,
                                                2.0 * Metre)), 4));
}

// These four tests cover all the branches of ToQuaternion.
TEST_F(RotationTest, ToQuaternion1) {
  R3Element<double> const v1 = {2, 5, 6};
  R3Element<double> const v2 =
      R3Element<double>({-3, 4, 1}).OrthogonalizationAgainst(v1);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix<double> m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion2) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> const v2 =
      R3Element<double>({-3, 4, 1}).OrthogonalizationAgainst(v1);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix<double> m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion3) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> const v2 =
      R3Element<double>({-3, 4, 1}).OrthogonalizationAgainst(v1);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix<double> m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationTest, ToQuaternion4) {
  R3Element<double> const v1 = {-2, -5, -6};
  R3Element<double> const v2 =
      R3Element<double>({-3, 4, 1}).OrthogonalizationAgainst(v1);
  R3Element<double> v3 = Cross(v1, v2);
  R3Element<double> const w1 = Normalize(v1);
  R3Element<double> const w2 = Normalize(v2);
  R3Element<double> const w3 = Normalize(v3);
  R3x3Matrix<double> m = {w1, w2, w3};
  Rot rotation(ToQuaternion(m.Transpose()));
  EXPECT_THAT(rotation(e1_).coordinates(), AlmostEquals(w1, 6));
  EXPECT_THAT(rotation(e2_).coordinates(), AlmostEquals(w2, 5));
  EXPECT_THAT(rotation(e3_).coordinates(), AlmostEquals(w3, 1));
}

TEST_F(RotationDeathTest, SerializationError) {
  Identity<World, World> id;
  EXPECT_DEATH({
    serialization::LinearMap message;
    id.WriteToMessage(&message);
    Rot const r = Rot::ReadFromMessage(message);
  }, "HasExtension.*Rotation");
}

TEST_F(RotationTest, SerializationSuccess) {
  serialization::LinearMap message;
  rotation_a_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_from_frame());
  EXPECT_TRUE(message.has_to_frame());
  EXPECT_EQ(message.from_frame().tag_type_fingerprint(),
            message.to_frame().tag_type_fingerprint());
  EXPECT_EQ(message.from_frame().tag(),
            message.to_frame().tag());
  EXPECT_EQ(message.from_frame().is_inertial(),
            message.to_frame().is_inertial());
  EXPECT_TRUE(message.HasExtension(serialization::Rotation::extension));
  serialization::Rotation const& extension =
      message.GetExtension(serialization::Rotation::extension);
  EXPECT_THAT(extension.quaternion().real_part(), AlmostEquals(0.5, 1));
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().x().double_());
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().y().double_());
  EXPECT_EQ(0.5, extension.quaternion().imaginary_part().z().double_());
  Rot const r = Rot::ReadFromMessage(message);
  EXPECT_EQ(rotation_a_(vector_), r(vector_));
}

TEST_F(RotationTest, Basis) {
  Vector<double, World> a = Normalize(Vector<double, World>({1, 1, -1}));
  Vector<double, World> b = Normalize(Vector<double, World>({1, 0, 1}));
  Bivector<double, World> c = Wedge(a, b);

  Rotation<World, World1> const to_world1(a, b, c);
  EXPECT_THAT(to_world1(a),
              Componentwise(AlmostEquals(1, 1),
                            AlmostEquals(0, 0),
                            VanishesBefore(1, 0)));
  EXPECT_THAT(to_world1(b),
              Componentwise(AlmostEquals(0, 0),
                            AlmostEquals(1, 0),
                            VanishesBefore(1, 0)));
  EXPECT_THAT(to_world1(c),
              Componentwise(VanishesBefore(1, 2),
                            AlmostEquals(0, 0),
                            AlmostEquals(1, 0)));

  Rotation<World1, World> const to_world(a, b, c);
  EXPECT_THAT(to_world(Vector<double, World1>({1, 0, 0})), AlmostEquals(a, 1));
  EXPECT_THAT(to_world(Vector<double, World1>({0, 1, 0})), AlmostEquals(b, 2));
  EXPECT_THAT(to_world(Bivector<double, World1>({0, 0, 1})),
              AlmostEquals(c, 4));
}

TEST_F(RotationTest, Enums) {
  Angle const α = 30 * Degree;
  Angle const β = 20 * Degree;
  Angle const γ = 100 * Degree;
  Rotation<World, World1> const zxz_euler(α, β, γ,
                                          EulerAngles::ZXZ,
                                          DefinesFrame<World1>{});

  // Checks that using the convention |axes| for Euler angles, conjugated by the
  // given |permutation|, and with appropriate sign changes, is equivalent to
  // the ZXZ convention.
  auto const check_euler_angles = [this, &α, &β, &γ, &zxz_euler](
      EulerAngles axes, auto permutation) {
    constexpr Handedness permuted_handedness =
        std::is_same_v<decltype(permutation), EvenPermutation>
            ? Handedness::Right
            : Handedness::Left;
    using Permuted =
        Frame<enum class PermutedTag, Inertial, permuted_handedness>;
    using PermutedRotated =
        Frame<enum class PermutedRotatedTag, Inertial, permuted_handedness>;
    Permutation<World, Permuted> const σ(permutation);
    Permutation<PermutedRotated, World1> const τ =
        (Permutation<Permuted, PermutedRotated>::Identity() * σ *
         Permutation<World1, World>::Identity()).Inverse();
    Rotation<Permuted, PermutedRotated> const euler(
        σ.Determinant() * α,
        σ.Determinant() * β,
        σ.Determinant() * γ,
        axes,
        DefinesFrame<PermutedRotated>{});
    EXPECT_THAT(τ(euler(σ(e1_))), Eq(zxz_euler(e1_)));
    EXPECT_THAT(τ(euler(σ(e2_))), Eq(zxz_euler(e2_)));
    EXPECT_THAT(τ(euler(σ(e3_))), Eq(zxz_euler(e3_)));
  };

  check_euler_angles(EulerAngles::ZXZ, EvenPermutation::XYZ);
  check_euler_angles(EulerAngles::XYX, EvenPermutation::ZXY);
  check_euler_angles(EulerAngles::YZY, EvenPermutation::YZX);

  check_euler_angles(EulerAngles::ZYZ, OddPermutation::YXZ);
  check_euler_angles(EulerAngles::XZX, OddPermutation::ZYX);
  check_euler_angles(EulerAngles::YXY, OddPermutation::XZY);

  Rotation<World, World1> const xyz_cardano(α, β, γ,
                                            CardanoAngles::XYZ,
                                            DefinesFrame<World1>{});

  // Checks that using the convention |axes| for Cardano angles, conjugated by
  // the given |permutation|, and with appropriate sign changes, is equivalent
  // to the XYZ convention.
  auto const check_cardano_angles = [this, &α, &β, &γ, &xyz_cardano](
                                        CardanoAngles axes, auto permutation) {
    constexpr Handedness permuted_handedness =
        std::is_same_v<decltype(permutation), EvenPermutation>
            ? Handedness::Right
            : Handedness::Left;
    using Permuted =
        Frame<enum class PermutedTag, Inertial, permuted_handedness>;
    using PermutedRotated =
        Frame<enum class PermutedRotatedTag, Inertial, permuted_handedness>;
    Permutation<World, Permuted> const σ(permutation);
    Permutation<PermutedRotated, World1> const τ =
        (Permutation<Permuted, PermutedRotated>::Identity() * σ *
         Permutation<World1, World>::Identity()).Inverse();
    Rotation<Permuted, PermutedRotated> const cardano(
        σ.Determinant() * α,
        σ.Determinant() * β,
        σ.Determinant() * γ,
        axes,
        DefinesFrame<PermutedRotated>{});
    EXPECT_THAT(τ(cardano(σ(e1_))), Eq(xyz_cardano(e1_)));
    EXPECT_THAT(τ(cardano(σ(e2_))), Eq(xyz_cardano(e2_)));
    EXPECT_THAT(τ(cardano(σ(e3_))), Eq(xyz_cardano(e3_)));
  };

  check_cardano_angles(CardanoAngles::XYZ, EvenPermutation::XYZ);
  check_cardano_angles(CardanoAngles::YZX, EvenPermutation::ZXY);
  check_cardano_angles(CardanoAngles::ZXY, EvenPermutation::YZX);

  check_cardano_angles(CardanoAngles::XZY, OddPermutation::XZY);
  check_cardano_angles(CardanoAngles::ZYX, OddPermutation::ZYX);
  check_cardano_angles(CardanoAngles::YXZ, OddPermutation::YXZ);
}

TEST_F(RotationTest, EulerAngles) {
  // Angles defining an orbit.
  Angle const Ω = 30 * Degree;
  Angle const i = 20 * Degree;
  Angle const ω = 100 * Degree;

  // The frame in which the above elements are given.
  using Reference = Frame<enum class ReferenceTag>;
  // |Nodes| shares its z axis with |Reference|, and has the ascending node of
  // the orbit as its positive x direction.
  using Nodes = Frame<enum class NodesTag>;
  // |Plane| also has the ascending node of the orbit as its positive x
  // direction, and has the orbital plane as its xy plane (with z being the
  // positive orbit normal).
  using Plane = Frame<enum class PlaneTag>;
  // |Orbit| has its x axis towards the periapsis, and its z axis towards the
  // positive orbit normal.
  using Orbit = Frame<enum class OrbitTag>;

  Bivector<double, Reference> const celestial_pole({0, 0, 1});
  Bivector<double, Nodes> const ascending_node({1, 0, 0});
  Bivector<double, Plane> const orbit_normal({0, 0, 1});

  Rotation<Reference, Nodes> const to_nodes(Ω,
                                            celestial_pole,
                                            DefinesFrame<Nodes>{});
  Rotation<Nodes, Plane> const to_plane(i,
                                        ascending_node,
                                        DefinesFrame<Plane>{});
  Rotation<Plane, Orbit> const to_orbit(ω,
                                        orbit_normal,
                                        DefinesFrame<Orbit>{});

  Rotation<Reference, Orbit> const to_orbit_direct(Ω, i, ω,
                                                   EulerAngles::ZXZ,
                                                   DefinesFrame<Orbit>{});

  EXPECT_THAT(to_orbit_direct.quaternion(),
              AlmostEquals((to_orbit * to_plane * to_nodes).quaternion(), 2));
}

TEST_F(RotationTest, CardanoAngles) {
  using Ground = Frame<enum class GroundTag>;
  Vector<double, Ground> north({1, 0, 0});
  Vector<double, Ground> east({0, 1, 0});
  Vector<double, Ground> down({0, 0, 1});
  auto const up = -down;

  using Aircraft = Frame<enum class AircraftTag>;
  Vector<double, Aircraft> forward({1, 0, 0});
  Vector<double, Aircraft> right({0, 1, 0});
  Vector<double, Aircraft> bottom({0, 0, 1});

  Angle const heading = 1 * Degree;
  Angle const pitch = 5 * Degree;

  // Level flight North.
  Velocity<Ground> const v = north * 100 * Metre / Second;

  {
    // No roll.
    Angle const roll = 0 * Degree;
    Rotation<Aircraft, Ground> const to_ground(heading,
                                               pitch,
                                               roll,
                                               CardanoAngles::ZYX,
                                               DefinesFrame<Aircraft>{});
    Rotation<Ground, Aircraft> const to_aircraft(heading,
                                                 pitch,
                                                 roll,
                                                 CardanoAngles::ZYX,
                                                 DefinesFrame<Aircraft>{});

    // Positive pitch is up.
    EXPECT_THAT(InnerProduct(to_ground(forward), up), Gt(0));
    // Small heading is slightly East from North.
    EXPECT_THAT(InnerProduct(to_ground(forward), east), Gt(0));
    // No roll, the wings point to the horizon.
    EXPECT_THAT(InnerProduct(to_ground(right), down), Eq(0));

    Velocity<Aircraft> const v_aircraft = to_aircraft(v);
    auto const spherical_aircraft_velocity =
        v_aircraft.coordinates().ToSpherical();
    Angle const angle_of_attack = spherical_aircraft_velocity.latitude;
    Angle const sideslip = spherical_aircraft_velocity.longitude;

    EXPECT_THAT(angle_of_attack, Gt(0 * Degree));
    // Positive angle of attack results in positive z velocity.
    EXPECT_THAT(v_aircraft.coordinates().z, Gt(0 * Metre / Second));

    EXPECT_THAT(sideslip, Lt(0 * Degree));
  }

  {
    // The angle between forward and north.  This is almost certainly a terrible
    // formula for that, cf. "A Case Study of Bits Lost in Space", in Kahan's
    // "How Futile are Mindless Assessments of Roundoff in Floating-Point
    // Computation?".
    Angle const angle_to_north = ArcCos(Cos(heading) * Cos(pitch));
    // This eliminates sideslip.
    Angle const roll = ArcSin(Sin(heading) / Sin(angle_to_north));
    EXPECT_THAT(roll, Gt(0 * Degree));

    Rotation<Aircraft, Ground> const to_ground(heading,
                                               pitch,
                                               roll,
                                               CardanoAngles::ZYX,
                                               DefinesFrame<Aircraft>{});
    Rotation<Ground, Aircraft> const to_aircraft(heading,
                                                 pitch,
                                                 roll,
                                                 CardanoAngles::ZYX,
                                                 DefinesFrame<Aircraft>{});

    EXPECT_THAT(AngleBetween(to_ground(forward), north),
                AlmostEquals(angle_to_north, 52));

    // Positive pitch is up.
    EXPECT_THAT(InnerProduct(to_ground(forward), up), Gt(0));
    // Small heading is slightly East from North.
    EXPECT_THAT(InnerProduct(to_ground(forward), east), Gt(0));
    // Positive roll, the right wing points down.
    EXPECT_THAT(InnerProduct(to_ground(right), down), Gt(0));

    Velocity<Aircraft> const v_aircraft = to_aircraft(v);
    auto const spherical_aircraft_velocity =
        v_aircraft.coordinates().ToSpherical();
    Angle const angle_of_attack = spherical_aircraft_velocity.latitude;
    Angle const sideslip = spherical_aircraft_velocity.longitude;

    EXPECT_THAT(angle_of_attack, AlmostEquals(angle_to_north, 52));
    EXPECT_THAT(sideslip, VanishesBefore(1 * Radian, 1));
  }

  {
    // Flying sideways (positive roll, right wing points down).
    Angle const roll = π / 2 * Radian;
    Rotation<Aircraft, Ground> const to_ground(heading,
                                               pitch,
                                               roll,
                                               CardanoAngles::ZYX,
                                               DefinesFrame<Aircraft>{});
    Rotation<Ground, Aircraft> const to_aircraft(heading,
                                                 pitch,
                                                 roll,
                                                 CardanoAngles::ZYX,
                                                 DefinesFrame<Aircraft>{});

    // Positive pitch is up.
    EXPECT_THAT(InnerProduct(to_ground(forward), up), Gt(0));
    // Small heading is slightly East from North.
    EXPECT_THAT(InnerProduct(to_ground(forward), east), Gt(0));
    // Positive roll makes the right wing point down.
    EXPECT_THAT(AngleBetween(to_ground(right), down), AlmostEquals(pitch, 1));

    Velocity<Aircraft> const v_aircraft = to_aircraft(v);
    auto const spherical_aircraft_velocity =
        v_aircraft.coordinates().ToSpherical();
    Angle const angle_of_attack = spherical_aircraft_velocity.latitude;
    Angle const sideslip = spherical_aircraft_velocity.longitude;

    EXPECT_THAT(angle_of_attack, AlmostEquals(heading, 3, 5));
    EXPECT_THAT(sideslip, AlmostEquals(pitch, 0, 2));
  }
}

}  // namespace internal_rotation
}  // namespace geometry
}  // namespace principia
