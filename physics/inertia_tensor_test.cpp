
#include "physics/inertia_tensor.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {

using geometry::Barycentre;
using geometry::Displacement;
using geometry::Frame;
using geometry::Position;
using geometry::R3x3Matrix;
using quantities::Density;
using quantities::Length;
using quantities::Mass;
using quantities::Pow;
using quantities::si::Kilogram;
using quantities::si::Metre;

TEST(InertiaTensorTest, Abdulghany) {
  using F1 = Frame<serialization::Frame::TestTag,
                   serialization::Frame::TEST1,
                   /*frame_is_inertial=*/true>;
  using F2 = Frame<serialization::Frame::TestTag,
                   serialization::Frame::TEST2,
                   /*frame_is_inertial=*/true>;
  using F3 = Frame<serialization::Frame::TestTag,
                   serialization::Frame::TEST2,
                   /*frame_is_inertial=*/true>;

  Density const ρ = 13593 * Kilogram / Pow<3>(Metre);
  Length const r = 3 * Metre;
  Mass const cuboid_mass = 32 * ρ * Pow<3>(r);
  Mass const cylinder_mass = 8 * π * ρ * Pow<3>(r);

  Position<F1> const cuboid_centre_of_mass = F1::origin;
  Position<F1> const cylinder_centre_of_mass =
      F1::origin + Displacement<F1>({0 * r, 3 * r, 5 * r});
  Position<F1> const overall_centre_of_mass =
      Barycentre<Position<F1>, Mass>(
          {cuboid_centre_of_mass, cylinder_centre_of_mass},
          {cuboid_mass, cylinder_mass});

  InertiaTensor<F1> const inertia_tensor_cuboid(
      cuboid_mass,
      32 * ρ * Pow<5>(r) *
          R3x3Matrix<double>{{68.0 / 12.0, 0.0, 0.0},
                             {0.0, 8.0 / 12.0, 0.0},
                             {0.0, 0.0, 68.0 / 12.0}},
      cuboid_centre_of_mass);
  InertiaTensor<F2> const inertia_tensor_cylinder(
      cylinder_mass,
      8 * π * ρ * Pow<5>(r) *
          R3x3Matrix<double>{{67.0 / 12.0, 0.0, 0.0},
                             {0.0, 67.0 / 12.0, 0.0},
                             {0.0, 0.0, 1.0 / 2.0}},
      F2::origin);

  InertiaTensor<F3> itcu =
      inertia_tensor_cuboid.Translate<F3>(overall_centre_of_mass);
  //InertiaTensor<F3> itcy =
  //    inertia_tensor_cylinder.Translate<F3>(overall_centre_of_mass);
  //InertiaTensor<F3> ittot = itcu + itcy;
}

}  // namespace physics
}  // namespace principia
