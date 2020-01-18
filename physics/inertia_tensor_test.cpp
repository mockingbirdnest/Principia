
#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using geometry::Barycentre;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::Rotation;
using quantities::Density;
using quantities::Length;
using quantities::Mass;
using quantities::MomentOfInertia;
using quantities::Pow;
using quantities::SIUnit;
using quantities::si::Degree;
using quantities::si::Kilogram;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::VanishesBefore;
using testing_utilities::operator""_⑴;
using ::testing::Matcher;

// See https://en.wikipedia.org/wiki/List_of_moments_of_inertia.
class InertiaTensorTest : public ::testing::Test {
 protected:
  using Frame0 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST>;
  using Frame1 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST1>;
  using Frame2 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST2>;
  using Frame3 = Frame<serialization::Frame::TestTag,
                       Inertial,
                       Handedness::Right,
                       serialization::Frame::TEST3>;

  template<typename Frame>
  void CheckMomentsOfInertia(
      InertiaTensor<Frame> const& tensor,
      Matcher<R3Element<MomentOfInertia>> const& matcher) {
    using PrincipalAxesFrame =
        geometry::Frame<enum class PrincipalAxesFrameTag>;
    auto const principal_axes =
        tensor.template Diagonalize<PrincipalAxesFrame>();
    EXPECT_THAT(principal_axes.moments_of_inertia, matcher);
  }

  template<typename Frame>
  void CheckInertiaTensorCoordinate(InertiaTensor<Frame> const& tensor,
                                    int const row,
                                    int const column,
                                    Matcher<MomentOfInertia> const& matcher) {
    auto const form_coordinates = tensor.form_.coordinates();
    if (row == column) {
      int const row1 = (row + 2) % 3;
      int const column1 = (column + 2) % 3;
      int const row2 = (row + 4) % 3;
      int const column2 = (column + 4) % 3;
      EXPECT_THAT(
          form_coordinates(row1, column1) + form_coordinates(row2, column2),
          matcher) << row << " " << column;
    } else {
      EXPECT_THAT(-form_coordinates(row, column), matcher)
          << row << " " << column;
    }
  }
};

// A point mass at a certain distance of an axis.
TEST_F(InertiaTensorTest, PointMass) {
  Mass const mass = 3 * Kilogram;

  using CentreOfMass = Frame0;
  using GeneralPoint = Frame1;

  InertiaTensor<CentreOfMass> const inertia_tensor_centre_of_mass(
      mass,
      SIUnit<MomentOfInertia>() * R3x3Matrix<double>{{0, 0, 0},
                                                     {0, 0, 0},
                                                     {0, 0, 0}},
      CentreOfMass::origin);

  auto const displacement =
      Displacement<CentreOfMass>({1 * Metre, 2 * Metre, 4 * Metre});
  RigidTransformation<CentreOfMass, GeneralPoint> const translation(
      CentreOfMass::origin + displacement,
      GeneralPoint::origin,
      Identity<CentreOfMass, GeneralPoint>().Forget<OrthogonalMap>());
  InertiaTensor<GeneralPoint> const inertia_tensor_general_point =
      inertia_tensor_centre_of_mass.Transform(translation);

  // One of the principal axes goes through the centre of mass and the general
  // point, and it has no inertia.  The other principal axes are orthogonal and
  // they have the same inertia, which is the elementary inertia of a point with
  // respect to an axis.
  MomentOfInertia const moment_of_inertia = mass * displacement.Norm²();
  CheckMomentsOfInertia(
      inertia_tensor_general_point,
      Componentwise(VanishesBefore(moment_of_inertia, 2),
                    RelativeErrorFrom(moment_of_inertia, IsNear(2.9e-9_⑴)),
                    RelativeErrorFrom(moment_of_inertia, IsNear(2.9e-9_⑴))));
}

// A rod with respect to an axis going through its extremity
TEST_F(InertiaTensorTest, Rod) {
  Mass const mass = 3 * Kilogram;
  Length const length = 5 * Metre;

  using CentreOfMass = Frame0;
  using Extremity = Frame1;

  InertiaTensor<CentreOfMass> const inertia_tensor_centre_of_mass(
      mass,
      mass * Pow<2>(length) * R3x3Matrix<double>{{0, 0, 0},
                                                 {0, 1.0 / 12.0, 0},
                                                 {0, 0, 1.0 / 12.0}},
      CentreOfMass::origin);

  auto const displacement =
      Displacement<CentreOfMass>({length / 2.0, 0 * Metre, 0 * Metre});
  RigidTransformation<CentreOfMass, Extremity> const translation(
      CentreOfMass::origin + displacement,
      Extremity::origin,
      Identity<CentreOfMass, Extremity>().Forget<OrthogonalMap>());
  InertiaTensor<Extremity> const inertia_tensor_extremity =
      inertia_tensor_centre_of_mass.Transform(translation);

  // One of the principal axes goes through the axis of the rod, and it has no
  // inertia.  The other principal axes are orthogonal and they have the same
  // inertia.
  MomentOfInertia const moment_of_inertia = mass * Pow<2>(length) / 3.0;
  CheckMomentsOfInertia(
      inertia_tensor_extremity,
      Componentwise(VanishesBefore(moment_of_inertia, 1),
                    RelativeErrorFrom(moment_of_inertia, IsNear(4.1e-9_⑴)),
                    RelativeErrorFrom(moment_of_inertia, IsNear(4.1e-9_⑴))));
}

// This test is from A. R. Abdulghany, Generalization of parallel axis theorem
// for rotational inertia [Abd17].
TEST_F(InertiaTensorTest, Abdulghany) {
  constexpr MomentOfInertia zero;
  Density const ρ = 1/*3593*/ * Kilogram / Pow<3>(Metre);  // Hg.

  using CylinderCentreOfMass = Frame0;
  using CuboidCentreOfMassZ = Frame1;
  using CuboidCentreOfMassY = Frame2;
  using OverallCentreOfMass = Frame3;

  Length const cylinder_radius = 3 * Metre;
  Length const cylinder_height = 24 * Metre;
  Mass const cylinder_mass = π * ρ * Pow<2>(cylinder_radius) * cylinder_height;

  // The cylinder with its axis along z.
  MomentOfInertia const cylinder_axis_inertia =
      cylinder_mass * Pow<2>(cylinder_radius) / 2.0;
  MomentOfInertia const cylinder_orthogonal_inertia =
      cylinder_mass *
      (3.0 * Pow<2>(cylinder_radius) + Pow<2>(cylinder_height)) / 12.0;
  InertiaTensor<CylinderCentreOfMass> const
      cylinder_inertia_centre_of_mass(
          cylinder_mass,
          R3x3Matrix<MomentOfInertia>{{cylinder_orthogonal_inertia, zero, zero},
                                      {zero, cylinder_orthogonal_inertia, zero},
                                      {zero, zero, cylinder_axis_inertia}},
          CylinderCentreOfMass::origin);

  Length const cuboid_small_side = 6 * Metre;
  Length const cuboid_long_side = 24 * Metre;
  Mass const cuboid_mass = ρ * Pow<2>(cuboid_small_side) * cuboid_long_side;

  // The cuboid with its long axis along z.
  MomentOfInertia const cuboid_long_axis_inertia =
      cuboid_mass * (Pow<2>(cuboid_small_side) + Pow<2>(cuboid_small_side)) /
      12.0;
  MomentOfInertia const cuboid_short_axis_inertia =
      cuboid_mass * (Pow<2>(cuboid_small_side) + Pow<2>(cuboid_long_side)) /
      12.0;
  InertiaTensor<CuboidCentreOfMassZ> const cuboid_inertia_centre_of_mass_z(
      cuboid_mass,
      R3x3Matrix<MomentOfInertia>{{cuboid_short_axis_inertia, zero, zero},
                                  {zero, cuboid_short_axis_inertia, zero},
                                  {zero, zero, cuboid_long_axis_inertia}},
      CuboidCentreOfMassZ::origin);

  // Rotate the cuboid around the x axis.
  RigidTransformation<CuboidCentreOfMassZ, CuboidCentreOfMassY> const
      cuboid_rotation(CuboidCentreOfMassZ::origin,
                      CuboidCentreOfMassY::origin,
                      Rotation<CuboidCentreOfMassZ, CuboidCentreOfMassY>(
                          90 * Degree,
                          Bivector<double, CuboidCentreOfMassZ>({1, 0, 0}),
                          DefinesFrame<CuboidCentreOfMassY>{})
                          .Forget<OrthogonalMap>());
  InertiaTensor<CuboidCentreOfMassY> const cuboid_inertia_centre_of_mass_y =
      cuboid_inertia_centre_of_mass_z.Transform(cuboid_rotation);

  // Determine the location of the overall centre of mass in both frames.
  R3Element<Length> const cylinder_cuboid_displacement(
      0 * Metre,
      cuboid_long_side / 2.0 - cylinder_radius,
      cylinder_height / 2.0 + cuboid_small_side / 2.0);
  Position<CuboidCentreOfMassY> const overall_centre_of_mass_cuboid =
      Barycentre<Position<CuboidCentreOfMassY>, Mass>(
          {CuboidCentreOfMassY::origin,
           CuboidCentreOfMassY::origin +
               Displacement<CuboidCentreOfMassY>(cylinder_cuboid_displacement)},
          {cuboid_mass, cylinder_mass});
  Position<CylinderCentreOfMass> const overall_centre_of_mass_cylinder =
      Barycentre<Position<CylinderCentreOfMass>, Mass>(
          {CylinderCentreOfMass::origin -
               Displacement<CylinderCentreOfMass>(cylinder_cuboid_displacement),
           CylinderCentreOfMass::origin},
          {cuboid_mass, cylinder_mass});

  // Compute the inertia of both objects with respect to the overall centre of
  // mass.
  InertiaTensor<OverallCentreOfMass> const
  cuboid_inertia_overall_centre_of_mass =
      cuboid_inertia_centre_of_mass_y.Transform(
          RigidTransformation<CuboidCentreOfMassY, OverallCentreOfMass>(
              overall_centre_of_mass_cuboid,
              OverallCentreOfMass::origin,
              Identity<CuboidCentreOfMassY, OverallCentreOfMass>()
              .Forget<OrthogonalMap>()));
  InertiaTensor<OverallCentreOfMass> const
  cylinder_inertia_overall_centre_of_mass =
      cylinder_inertia_centre_of_mass.Transform(
          RigidTransformation<CylinderCentreOfMass, OverallCentreOfMass>(
              overall_centre_of_mass_cylinder,
              OverallCentreOfMass::origin,
              Identity<CylinderCentreOfMass, OverallCentreOfMass>()
              .Forget<OrthogonalMap>()));

  // Finally, the overall inertia at the overall centre of mass.
  InertiaTensor<OverallCentreOfMass> const
      overall_inertia_overall_centre_of_mass =
          cuboid_inertia_overall_centre_of_mass +
          cylinder_inertia_overall_centre_of_mass;

  // The expectations come from [Abd17].
  CheckInertiaTensorCoordinate(
      overall_inertia_overall_centre_of_mass,
      0,
      0,
      AlmostEquals((162.0 * (1088.0 + 2172.0 * π + 67.0 * Pow<2>(π))) /
                       (4.0 + π) * SIUnit<MomentOfInertia>(),
                   1));
  CheckInertiaTensorCoordinate(overall_inertia_overall_centre_of_mass,
                               0,
                               1,
                               AlmostEquals(0 * SIUnit<MomentOfInertia>(), 0));
  CheckInertiaTensorCoordinate(overall_inertia_overall_centre_of_mass,
                               0,
                               2,
                               AlmostEquals(0 * SIUnit<MomentOfInertia>(), 0));
  CheckInertiaTensorCoordinate(overall_inertia_overall_centre_of_mass,
                               1,
                               0,
                               AlmostEquals(0 * SIUnit<MomentOfInertia>(), 0));
  CheckInertiaTensorCoordinate(
      overall_inertia_overall_centre_of_mass,
      1,
      1,
      AlmostEquals((162.0 * (128.0 + 1500.0 * π + 67.0 * Pow<2>(π))) /
                       (4.0 + π) * SIUnit<MomentOfInertia>(),
                   2));
  CheckInertiaTensorCoordinate(
      overall_inertia_overall_centre_of_mass,
      1,
      2,
      AlmostEquals(-116640.0 * π / (4.0 + π) * SIUnit<MomentOfInertia>(),
                   0, 1));
  CheckInertiaTensorCoordinate(overall_inertia_overall_centre_of_mass,
                               2,
                               0,
                               AlmostEquals(0 * SIUnit<MomentOfInertia>(), 0));
  CheckInertiaTensorCoordinate(
      overall_inertia_overall_centre_of_mass,
      2,
      1,
      AlmostEquals(-116640.0 * π / (4.0 + π) * SIUnit<MomentOfInertia>(),
                   0, 1));
  CheckInertiaTensorCoordinate(
      overall_inertia_overall_centre_of_mass,
      2,
      2,
      AlmostEquals((324.0 * (544.0 + 364.0 * π + 3.0 * Pow<2>(π))) / (4.0 + π) *
                       SIUnit<MomentOfInertia>(),
                   0));
}

}  // namespace internal_inertia_tensor
}  // namespace physics
}  // namespace principia
