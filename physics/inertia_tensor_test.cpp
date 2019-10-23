
#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {

using geometry::Barycentre;
using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::R3Element;
using geometry::R3x3Matrix;
using quantities::Density;
using quantities::Length;
using quantities::Mass;
using quantities::MomentOfInertia;
using quantities::Pow;
using quantities::SIUnit;
using quantities::si::Kilogram;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;

// See https://en.wikipedia.org/wiki/List_of_moments_of_inertia.
class InertiaTensorTest : public ::testing::Test {
 protected:
  using Frame1 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST1,
                       /*frame_is_inertial=*/true>;
  using Frame2 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST2,
                       /*frame_is_inertial=*/true>;
  using Frame3 = Frame<serialization::Frame::TestTag,
                       serialization::Frame::TEST3,
                       /*frame_is_inertial=*/true>;

  template<typename Frame>
  void CheckMomentsOfInertia(
      InertiaTensor<Frame> const& tensor,
      R3Element<MomentOfInertia> const& moments_of_inertia) {
    struct PrincipalAxesFrame {};

    auto const principal_axes = tensor.Diagonalize<PrincipalAxesFrame>();
    EXPECT_THAT(principal_axes.moments_of_inertia,
                AlmostEquals(moments_of_inertia, 0));
  }
};

TEST_F(InertiaTensorTest, PointMass) {
  Mass const mass = 3 * Kilogram;

  using CentreOfMass = Frame1;
  using GeneralPoint = Frame2;

  InertiaTensor<CentreOfMass> const inertia_tensor_centre_of_mass(
      mass,
      SIUnit<MomentOfInertia>() * R3x3Matrix<double>{{0, 0, 0},
                                                     {0, 0, 0},
                                                     {0, 0, 0}},
      CentreOfMass::origin);
  InertiaTensor<GeneralPoint> const inertia_tensor_general_point =
      inertia_tensor_centre_of_mass.Translate<GeneralPoint>(
          CentreOfMass::origin +
          Displacement<CentreOfMass>({1 * Metre, 2 * Metre, 4 * Metre}));

  CheckMomentsOfInertia(
      inertia_tensor_general_point,
      SIUnit<MomentOfInertia>() * R3Element<double>(0, 0, 63));
}

//TEST(InertiaTensorTest, Abdulghany) {
//  using F1 = Frame<serialization::Frame::TestTag,
//                   serialization::Frame::TEST1,
//                   /*frame_is_inertial=*/true>;
//  using F2 = Frame<serialization::Frame::TestTag,
//                   serialization::Frame::TEST2,
//                   /*frame_is_inertial=*/true>;
//  using F3 = Frame<serialization::Frame::TestTag,
//                   serialization::Frame::TEST2,
//                   /*frame_is_inertial=*/true>;
//
//  Density const ρ = 13593 * Kilogram / Pow<3>(Metre);
//  Length const r = 3 * Metre;
//  Mass const cuboid_mass = 32 * ρ * Pow<3>(r);
//  Mass const cylinder_mass = 8 * π * ρ * Pow<3>(r);
//
//  Position<F1> const cuboid_centre_of_mass = F1::origin;
//  Position<F1> const cylinder_centre_of_mass =
//      F1::origin + Displacement<F1>({0 * r, 3 * r, 5 * r});
//  Position<F1> const overall_centre_of_mass =
//      Barycentre<Position<F1>, Mass>(
//          {cuboid_centre_of_mass, cylinder_centre_of_mass},
//          {cuboid_mass, cylinder_mass});
//
//  InertiaTensor<F1> const inertia_tensor_cuboid(
//      cuboid_mass,
//      32 * ρ * Pow<5>(r) *
//          R3x3Matrix<double>{{68.0 / 12.0, 0.0, 0.0},
//                             {0.0, 8.0 / 12.0, 0.0},
//                             {0.0, 0.0, 68.0 / 12.0}},
//      cuboid_centre_of_mass);
//  InertiaTensor<F2> const inertia_tensor_cylinder(
//      cylinder_mass,
//      8 * π * ρ * Pow<5>(r) *
//          R3x3Matrix<double>{{67.0 / 12.0, 0.0, 0.0},
//                             {0.0, 67.0 / 12.0, 0.0},
//                             {0.0, 0.0, 1.0 / 2.0}},
//      F2::origin);
//
//  InertiaTensor<F3> itcu =
//      inertia_tensor_cuboid.Translate<F3>(overall_centre_of_mass);
//  InertiaTensor<F3> itcy =
//      inertia_tensor_cylinder.Translate<F3>(overall_centre_of_mass);
//  InertiaTensor<F3> ittot = itcu + itcy;
//}

}  // namespace physics
}  // namespace principia
