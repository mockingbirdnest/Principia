#include "mathematica/mathematica.hpp"

#include <list>
#include <string>
#include <vector>

#include "astronomy/orbital_elements.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "gtest/gtest.h"
#include "numerics/fixed_arrays.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/unbounded_arrays.hpp"
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
using geometry::Interval;
using geometry::Point;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::SymmetricBilinearForm;
using geometry::Vector;
using geometry::Velocity;
using numerics::FixedVector;
using numerics::HornerEvaluator;
using numerics::PiecewisePoissonSeries;
using numerics::PoissonSeries;
using numerics::PolynomialInMonomialBasis;
using numerics::UnboundedLowerTriangularMatrix;
using numerics::UnboundedUpperTriangularMatrix;
using numerics::UnboundedVector;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Length;
using quantities::Speed;
using quantities::Time;
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
        "List[" +
            (ToMathematica(1 * Metre) + "," + ToMathematica(2 * Second) + "," +
             ToMathematica(3 * Metre / Second)) +
            "]",
        ToMathematica(std::tuple{1 * Metre, 2 * Second, 3 * Metre / Second}));
  }
  {
    EXPECT_EQ(ToMathematica(std::tuple{2.0, 3.0}),
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
    EXPECT_EQ(ToMathematica(std::tuple{2.0, 3.0}),
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
    EXPECT_EQ(ToMathematica(std::tuple{2.0, 3.0}),
              ToMathematica(FixedVector<double, 2>({2, 3})));
    EXPECT_EQ("List[]", ToMathematica(FixedVector<int, 0>()));
  }
  {
    EXPECT_EQ(ToMathematica(std::tuple{2.0, 3.0, -4.0}),
              ToMathematica(R3Element<double>(2.0, 3.0, -4.0)));
  }
  {
    EXPECT_EQ(ToMathematica(std::tuple{std::tuple{1.0, 2.0, 3.0},
                                       std::tuple{4.0, 5.0, 6.0},
                                       std::tuple{7.0, 8.0, 9.0}}),
              ToMathematica(R3x3Matrix<double>({1.0, 2.0, 3.0},
                                               {4.0, 5.0, 6.0},
                                               {7.0, 8.0, 9.0})));
  }
  {
    EXPECT_EQ(
        "Quaternion[" +
            (ToMathematica(1.0) + "," + ToMathematica(2.0) + "," +
             ToMathematica(3.0) + "," + ToMathematica(-4.0)) +
            "]",
        ToMathematica(Quaternion(1.0, R3Element<double>(2.0, 3.0, -4.0))));
  }
  {
    EXPECT_EQ("Quantity[" + ToMathematica(1.5) + ",\" m s^-1\"]",
              ToMathematica(1.5 * Metre / Second));
  }
  {
    Vector<double, F> const v({2.0, 3.0, -4.0});
    EXPECT_EQ(ToMathematica(v.coordinates()), ToMathematica(v));
  }
  {
    Bivector<double, F> const b({2.0, 3.0, -4.0});
    EXPECT_EQ(ToMathematica(b.coordinates()), ToMathematica(b));
  }
  {
    Vector<double, F> const v({2.0, 3.0, -4.0});
    EXPECT_EQ(ToMathematica(v),
              ToMathematica(Point<Vector<double, F>>() + v));
  }
  {
    EXPECT_EQ(
        ToMathematica(SymmetricBilinearForm<double, F, Vector>().coordinates()),
        ToMathematica(SymmetricBilinearForm<double, F, Vector>()));
  }
  {
    DegreesOfFreedom<F> const dof(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second}));
    EXPECT_EQ(ToMathematica(std::tuple{dof.position(), dof.velocity()}),
              ToMathematica(dof));
  }
  {
    UnboundedLowerTriangularMatrix<double> l2({1,
                                               2, 3});
    EXPECT_EQ(
        "List[List[" + ToMathematica(1.0) + "," + ToMathematica(0.0) +
        "],List[" + ToMathematica(2.0) + "," + ToMathematica(3.0) + "]]",
        ToMathematica(l2));
  }
  {
    UnboundedUpperTriangularMatrix<double> u2({1, 2,
                                                  3});
    EXPECT_EQ(ToMathematica(std::tuple{std::tuple{1.0, 2.0},
                                       std::tuple{0.0, 3.0}}),
              ToMathematica(u2));
  }
  {
    UnboundedVector<double> v2({1, 2});
    EXPECT_EQ(ToMathematica(std::tuple{1.0, 2.0}), ToMathematica(v2));
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
    EXPECT_EQ(ToMathematica(std::tuple{trajectory.front().time,
                                       trajectory.front().degrees_of_freedom}),
              ToMathematica(trajectory.front()));
  }
  {
    OrbitalElements::EquinoctialElements elements{
        Instant(), 1 * Metre, 2, 3, 4 * Radian, 5, 6, 7, 8};
    EXPECT_EQ(
        ToMathematica(std::tuple{
            Instant(), 1 * Metre, 2.0, 3.0, 4 * Radian, 5.0, 6.0, 7.0, 8.0}),
        ToMathematica(elements));
  }
  {
    PolynomialInMonomialBasis<Length, Time, 2, HornerEvaluator> polynomial1(
        {2 * Metre, -3 * Metre / Second, 4 * Metre / Second / Second});
    EXPECT_EQ("Function[Plus[" + ToMathematica(2 * Metre) + ",Times[" +
                  ToMathematica(-3 * Metre / Second) + ",#],Times[" +
                  ToMathematica(4 * Metre / Second / Second) + ",Power[#,2]]]]",
        ToMathematica(polynomial1));
    PolynomialInMonomialBasis<Length, Instant, 2, HornerEvaluator> polynomial2(
        {5 * Metre, 6 * Metre / Second, -7 * Metre / Second / Second},
        Instant() + 2 * Second);
    EXPECT_EQ("Function[Plus[" + ToMathematica(5 * Metre) + ",Times[" +
                  ToMathematica(6 * Metre / Second) + ",Subtract[#," +
                  ToMathematica(polynomial2.origin()) + "]],Times[" +
                  ToMathematica(-7 * Metre / Second / Second) +
                  ",Power[Subtract[#," + ToMathematica(polynomial2.origin()) +
                  "],2]]]]",
              ToMathematica(polynomial2));
  }
  {
    using Series = PoissonSeries<double, 0, 0, HornerEvaluator>;
    Instant const t0 = Instant() + 3 * Second;
    Series::AperiodicPolynomial secular({1.5}, t0);
    Series::PeriodicPolynomial sin({0.5}, t0);
    Series::PeriodicPolynomial cos({-1}, t0);
    Series series(secular, {{4 * Radian / Second, {sin, cos}}});
    EXPECT_EQ(
        "Function[" +
            ("Plus[" + ToMathematicaBody(secular) + "," +
             ("Times[" + ToMathematicaBody(sin) + "," +
              ("Sin[Times[" + ToMathematica(4 * Radian / Second) + "," +
               ("Subtract[#," + ToMathematica(series.origin()) + "]") + "]]") +
              "]") + "," +
             ("Times[" + ToMathematicaBody(cos) + "," +
              ("Cos[Times[" + ToMathematica(4 * Radian / Second) + "," +
               ("Subtract[#," + ToMathematica(series.origin()) + "]") + "]]") +
              "]") +
             "]") +
            "]",
        ToMathematica(series));
  }
  {
    using PiecewiseSeries =
        PiecewisePoissonSeries<double, 0, 0, HornerEvaluator>;
    using Series = PiecewiseSeries::Series;
    Instant const t0 = Instant() + 3 * Second;
    Series series(Series::AperiodicPolynomial({1.5}, t0),
                  {{4 * Radian / Second,
                    {/*sin=*/Series::PeriodicPolynomial({0.5}, t0),
                     /*cos=*/Series::PeriodicPolynomial({-1}, t0)}}});
    Interval<Instant> interval{t0 - 2 * Second, t0 + 3 * Second};
    PiecewiseSeries pw(interval, series);
    EXPECT_EQ(
        "Function[Piecewise[List[List[" +
            (ToMathematicaBody(series) + "," +
             ("Between[#," +
              ToMathematica(std::tuple{interval.min, interval.max}) + "]")) +
            "]]]]",
        ToMathematica(pw));
  }
  {
    std::optional<std::string> opt1;
    std::optional<std::string> opt2("foo");
    EXPECT_EQ("List[]", ToMathematica(opt1));
    EXPECT_EQ("List[\"foo\"]", ToMathematica(opt2));
  }
}

TEST_F(MathematicaTest, Option) {
  EXPECT_EQ("Rule[option," + ToMathematica(3.0) + "]", Option("option", 3.0));
}

TEST_F(MathematicaTest, Assign) {
  EXPECT_EQ("Set[var," + ToMathematica(3.0) + "];\n", Assign("var", 3.0));
}

TEST_F(MathematicaTest, PlottableDataset) {
  std::vector<double> const abscissæ{2, 3};
  std::vector<std::string> const ordinates{"2", "3"};
  EXPECT_EQ("Transpose[" + ToMathematica(std::tuple{abscissæ, ordinates}) + "]",
            PlottableDataset(abscissæ, ordinates));
}

TEST_F(MathematicaTest, Escape) {
  EXPECT_EQ(R"("foo")", ToMathematica("foo"));
  {
    // This string messes up the macro.
    std::string expected = R"("fo\"o")";
    EXPECT_EQ(expected, ToMathematica("fo\"o"));
  }
  EXPECT_EQ(R"("fo\\o")", ToMathematica("fo\\o"));
  EXPECT_EQ(R"("")", ToMathematica(""));
}

TEST_F(MathematicaTest, ExpressIn) {
  {
    EXPECT_EQ(
        ToMathematica(3.0),
        ToMathematica(3.0 * Metre / Second / Second, ExpressIn(Metre, Second)));
    EXPECT_EQ(ToMathematica(1 * Radian / (1 * Degree)),
              ToMathematica(1 * Radian, ExpressIn(Degree)));
  }
  {
    Vector<Speed, F> const v(
        {2.0 * Metre / Second, 3.0 * Metre / Second, -4.0 * Metre / Second});
    EXPECT_EQ(ToMathematica(v / (Metre / Second)),
              ToMathematica(v, ExpressIn(Metre, Second)));
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
        ToMathematica(std::tuple{0.0,
                                 std::tuple{std::tuple{2.0, 3.0, -4.0},
                                            std::tuple{-1.0, -5.0, 8.0}}}),
        ToMathematica(trajectory.front(), ExpressIn(Metre, Second)));
  }
  {
    OrbitalElements::EquinoctialElements elements{
        Instant(), 1 * Metre, 2, 3, 4 * Radian, 5, 6, 7, 8};
    EXPECT_EQ(
        ToMathematica(std::tuple(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)),
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
    Logger logger(TEMP_DIR / "mathematica_test.wl");
    logger.Append("a", std::vector{1.0, 2.0, 3.0});
    logger.Append(u8"β", 4 * Metre / Second);
    logger.Append("a", F::origin);
    logger.Set("c", 5.0);
  }
  // Go check the file.
  EXPECT_EQ(Assign("a", std::tuple{std::vector{1.0, 2.0, 3.0}, F::origin}) +
                Assign(u8"β", std::tuple{4 * Metre / Second}) +
                Assign("c", 5.0),
            (std::stringstream{}
             << std::ifstream(TEMP_DIR / "mathematica_test0.wl").rdbuf())
                .str());
}
#endif

}  // namespace mathematica
}  // namespace principia
