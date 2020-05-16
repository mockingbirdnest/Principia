#include "mathematica/mathematica.hpp"

#include <list>
#include <string>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {

using astronomy::OrbitalElements;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::Point;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::Vector;
using geometry::Velocity;
using numerics::FixedVector;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Speed;
using quantities::si::Degree;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

class MathematicaTest : public ::testing::Test {
 protected:
  using F = Frame<enum class FTag>;
};

TEST_F(MathematicaTest, ToMathematica) {
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(std::vector<double>{2, 3}));
    EXPECT_EQ("List[]", ToMathematica(std::vector<int>{}));
  }
  {
    std::list<double> list;
    list.push_back(3.0);
    list.push_front(2.0);
    list.push_back(4.0);
    auto it = list.cbegin();
    ++it;
    ++it;
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(list.cbegin(), it));
    EXPECT_EQ("List[]", ToMathematica(it, it));
  }
  {
    EXPECT_EQ("True", ToMathematica(true));
    EXPECT_EQ("False", ToMathematica(false));
  }
  {
    EXPECT_EQ("3", ToMathematica(3));
    EXPECT_EQ("-2", ToMathematica(-2));
  }
  {
    EXPECT_EQ("SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]",
              ToMathematica(3.0));
    EXPECT_EQ("SetPrecision[-2.00000000000000000*^+09,$MachinePrecision]",
              ToMathematica(-2.0e9));
    EXPECT_EQ("SetPrecision[-0.00000000000000000*^+00,$MachinePrecision]",
              ToMathematica(-0.0));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(FixedVector<double, 2>({2, 3})));
    EXPECT_EQ("List[]", ToMathematica(FixedVector<int, 0>()));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(R3Element<double>(2.0, 3.0, -4.0)));
  }
  {
    EXPECT_EQ(
        "List["
        "List["
        "SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]],"
        "List["
        "SetPrecision[+4.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+5.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+6.00000000000000000*^+00,$MachinePrecision]],"
        "List["
        "SetPrecision[+7.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+9.00000000000000000*^+00,$MachinePrecision]]]",
        ToMathematica(R3x3Matrix<double>({1.0, 2.0, 3.0},
                                         {4.0, 5.0, 6.0},
                                         {7.0, 8.0, 9.0})));
  }
  {
    EXPECT_EQ(
        "Quaternion["
        "SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(Quaternion(1.0, R3Element<double>(2.0, 3.0, -4.0))));
  }
  {
    EXPECT_EQ(
        "Quantity["
        "SetPrecision[+1.50000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"]",
        ToMathematica(1.5 * Metre / Second));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(Vector<double, F>({2.0, 3.0, -4.0})));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(Bivector<double, F>({2.0, 3.0, -4.0})));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(Point<Vector<double, F>>() +
                      Vector<double, F>({2.0, 3.0, -4.0})));
  }
  {
    EXPECT_EQ(
        "List["
        "List["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision]],"
        "List["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision]],"
        "List["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision]]]",
        ToMathematica(SymmetricBilinearForm<double, F, Vector>()));
  }
  {
    EXPECT_EQ(
        "List["
        "List["
        "Quantity["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "Quantity["
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "Quantity["
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"]],"
        "List["
        "Quantity["
        "SetPrecision[-1.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"],"
        "Quantity["
        "SetPrecision[-5.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"],"
        "Quantity["
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"]]]",
        ToMathematica(DegreesOfFreedom<F>(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second}))));
  }
  {
    EXPECT_EQ(
        "List["
        "Quantity["
        "SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "Quantity["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "\" s\"],"
        "Quantity["
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"]]",
        ToMathematica(std::tuple{1 * Metre, 2 * Second, 3 * Metre / Second}));
  }
  {
    DiscreteTrajectory<F> trajectory;
    trajectory.Append(
        Instant(),
        DegreesOfFreedom<F>(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second})));
    EXPECT_EQ(
        "List["
        "Quantity["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "\" s\"],"
        "List["
        "List["
        "Quantity["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "Quantity["
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "Quantity["
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"]],"
        "List["
        "Quantity["
        "SetPrecision[-1.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"],"
        "Quantity["
        "SetPrecision[-5.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"],"
        "Quantity["
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision],"
        "\" m s^-1\"]]]]",
        ToMathematica(*trajectory.begin()));
  }
  {
    OrbitalElements::EquinoctialElements elements{
        Instant(), 1 * Metre, 2, 3, 4 * Radian, 5, 6, 7, 8};
    EXPECT_EQ(
        "List["
        "Quantity["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "\" s\"],"
        "Quantity[SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
        "\" m\"],"
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "Quantity[SetPrecision[+4.00000000000000000*^+00,$MachinePrecision],"
        "\" rad\"],"
        "SetPrecision[+5.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+6.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+7.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(elements));
  }
  {
    std::optional<std::string> opt1;
    std::optional<std::string> opt2("foo");
    EXPECT_EQ("List[]", ToMathematica(opt1));
    EXPECT_EQ("List[foo]", ToMathematica(opt2));
  }
  {
    EXPECT_EQ("foo\"bar", ToMathematica("foo\"bar"));
  }
}

TEST_F(MathematicaTest, Apply) {
  EXPECT_EQ("Fun[]", Apply("Fun", {}));
  EXPECT_EQ("Fun[1,[],\"string\"]", Apply("Fun", {"1", "[]", "\"string\""}));
}

TEST_F(MathematicaTest, Option) {
  EXPECT_EQ(
      "Rule["
      "option,"
      "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]]",
      Option("option", ToMathematica(3.0)));
}

TEST_F(MathematicaTest, Assign) {
  EXPECT_EQ(
      "Set["
      "var,"
      "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]];\n",
      Assign("var", ToMathematica(3.0)));
}

TEST_F(MathematicaTest, PlottableDataset) {
  EXPECT_EQ(
      "Transpose["
      "List["
      "List["
      "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
      "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]],"
      "List[2,3]]]",
      PlottableDataset(std::vector<double>{2, 3},
                       std::vector<std::string>{"2", "3"}));
}

TEST_F(MathematicaTest, Escape) {
  EXPECT_EQ(R"("foo")", Escape("foo"));
  {
    // This string messes up the macro.
    std::string expected = R"("fo\"o")";
    EXPECT_EQ(expected, Escape("fo\"o"));
  }
  EXPECT_EQ(R"("fo\\o")", Escape("fo\\o"));
  EXPECT_EQ(R"("")", Escape(""));
}

TEST_F(MathematicaTest, ExpressIn) {
  {
    EXPECT_EQ("SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]",
              ToMathematica(3.0 * Metre / Second / Second,
                            ExpressIn(Metre, Second)));
    EXPECT_EQ("SetPrecision[+5.72957795130823229*^+01,$MachinePrecision]",
              ToMathematica(1 * Radian, ExpressIn(Degree)));
  }
  {
    EXPECT_EQ(
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(Vector<Speed, F>({2.0 * Metre / Second,
                                        3.0 * Metre / Second,
                                        -4.0 * Metre / Second}),
                      ExpressIn(Metre, Second)));
  }
  {
    DiscreteTrajectory<F> trajectory;
    trajectory.Append(
        Instant(),
        DegreesOfFreedom<F>(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second})));
    EXPECT_EQ(
        "List["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "List["
        "List["
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-4.00000000000000000*^+00,$MachinePrecision]],"
        "List["
        "SetPrecision[-1.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[-5.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision]]]]",
        ToMathematica(*trajectory.begin(), ExpressIn(Metre, Second)));
  }
  {
    OrbitalElements::EquinoctialElements elements{
        Instant(), 1 * Metre, 2, 3, 4 * Radian, 5, 6, 7, 8};
    EXPECT_EQ(
        "List["
        "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+4.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+5.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+6.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+7.00000000000000000*^+00,$MachinePrecision],"
        "SetPrecision[+8.00000000000000000*^+00,$MachinePrecision]]",
        ToMathematica(elements, ExpressIn(Metre, Second, Radian)));
  }

// Does not compile, by design.
#if 0
  ToMathematica(1 * Radian, "foobar");
#endif
#if 0
  ToMathematica(1 * Radian, ExpressIn(Metre));
#endif
#if 0
  ToMathematica(1 * Radian, ExpressIn(Degree, Metre, Metre));
#endif
}

// std::filesystem is broken on macOS.
#if !defined(__APPLE__)
TEST_F(MathematicaTest, Logger) {
  {
    Logger logger(TEMP_DIR / "mathematica_test.wl", /*make_unique=*/true);
    logger.Append("a", std::vector{1.0, 2.0, 3.0});
    logger.Append("b", 4 * Metre / Second);
    logger.Append("a", F::origin);
  }
  // Go check the file.
  EXPECT_EQ(
      "Set["
      "a,"
      "List["
      "List["
      "SetPrecision[+1.00000000000000000*^+00,$MachinePrecision],"
      "SetPrecision[+2.00000000000000000*^+00,$MachinePrecision],"
      "SetPrecision[+3.00000000000000000*^+00,$MachinePrecision]],"
      "List["
      "Quantity["
      "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
      "\" m\"],"
      "Quantity["
      "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
      "\" m\"],"
      "Quantity["
      "SetPrecision[+0.00000000000000000*^+00,$MachinePrecision],"
      "\" m\"]]]];\n"
      "Set["
      "b,"
      "List["
      "Quantity["
      "SetPrecision[+4.00000000000000000*^+00,$MachinePrecision],"
      "\" m s^-1\"]]];\n",
      (
          std::stringstream{}
          << std::ifstream(TEMP_DIR / "mathematica_test0.wl").rdbuf())
          .str());
}
#endif

}  // namespace mathematica
}  // namespace principia
