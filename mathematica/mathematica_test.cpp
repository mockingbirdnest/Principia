#include "mathematica/mathematica.hpp"

#include <list>
#include <string>
#include <vector>

#include "absl/strings/str_replace.h"
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
#include "physics/discrete_traject0ry.hpp"
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
using numerics::DoublePrecision;
using numerics::FixedVector;
using numerics::HornerEvaluator;
using numerics::PiecewisePoissonSeries;
using numerics::PoissonSeries;
using numerics::PolynomialInMonomialBasis;
using numerics::UnboundedLowerTriangularMatrix;
using numerics::UnboundedUpperTriangularMatrix;
using numerics::UnboundedVector;
using physics::DegreesOfFreedom;
using physics::DiscreteTraject0ry;
using quantities::Infinity;
using quantities::Length;
using quantities::Speed;
using quantities::Sqrt;
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
        absl::StrReplaceAll(u8"List[α,β,γ]",
                            {{u8"α", ToMathematica(1 * Metre)},
                             {u8"β", ToMathematica(2 * Second)},
                             {u8"γ", ToMathematica(3 * Metre / Second)}}),
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
    EXPECT_EQ("Times[16^^13456789ABCDEF,Power[2,Subtract[-163,52]]]",
              ToMathematica(0x1.3'4567'89AB'CDEFp-163));
    EXPECT_EQ("Minus[Times[16^^10000000000000,Power[2,Subtract[-1074,52]]]]",
              ToMathematica(-0x1p-1074));
    EXPECT_EQ("Minus[Times[16^^1FFFFFFFFFFFFF,Power[2,Subtract[1023,52]]]]",
              ToMathematica(-0x1.F'FFFF'FFFF'FFFFp1023));
    EXPECT_EQ("Times[16^^19ABCDE,Power[2,Subtract[-12,24]]]",
              ToMathematica(0x1.9ABCDEp-12f));
    EXPECT_EQ("0", ToMathematica(0.0));
    EXPECT_EQ("Minus[0]", ToMathematica(-0.0));  // Not that this does anything.
    EXPECT_EQ("Infinity", ToMathematica(Infinity<double>));
    EXPECT_EQ("Minus[Infinity]", ToMathematica(-Infinity<double>));
    EXPECT_EQ("Minus[Indeterminate]", ToMathematica(Sqrt(-1)));
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
        absl::StrReplaceAll(u8"Quaternion[α,β,γ,δ]",
                            {{u8"α", ToMathematica(1.0)},
                             {u8"β", ToMathematica(2.0)},
                             {u8"γ", ToMathematica(3.0)},
                             {u8"δ", ToMathematica(-4.0)}}),
        ToMathematica(Quaternion(1.0, R3Element<double>(2.0, 3.0, -4.0))));
  }
  {
    EXPECT_EQ(absl::StrReplaceAll(u8R"(Quantity[α," m s^-1"])",
                                  {{u8"α", ToMathematica(1.5)}}),
              ToMathematica(1.5 * Metre / Second));
  }
  {
    DoublePrecision<double> d(3);
    d += 5e-20;
    EXPECT_EQ(absl::StrReplaceAll(
                  u8"Plus[α,β]",
                  {{u8"α", ToMathematica(3.0)},
                   {u8"β", ToMathematica(5e-20)}}),
              ToMathematica(d));
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
        ToMathematica(std::tuple{std::tuple{1.0, 0.0}, std::tuple{2.0, 3.0}}),
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
    DiscreteTraject0ry<F> trajectory;
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
    EXPECT_EQ(absl::StrReplaceAll(
                  u8"Function[Plus[α,Times[β,#],Times[γ,Power[#,2]]]]",
                  {{u8"α", ToMathematica(2 * Metre)},
                   {u8"β", ToMathematica(-3 * Metre / Second)},
                   {u8"γ", ToMathematica(4 * Metre / Second / Second)}}),
              ToMathematica(polynomial1));
    PolynomialInMonomialBasis<Length, Instant, 2, HornerEvaluator> polynomial2(
        {5 * Metre, 6 * Metre / Second, -7 * Metre / Second / Second},
        Instant() + 2 * Second);
    EXPECT_EQ(absl::StrReplaceAll(
                  u8R"(Function[Plus[
                        α,
                        Times[β,Subtract[#,δ]],
                        Times[γ,Power[Subtract[#,δ],2]]]])",
                  {{u8"α", ToMathematica(5 * Metre)},
                   {u8"β", ToMathematica(6 * Metre / Second)},
                   {u8"γ", ToMathematica(-7 * Metre / Second / Second)},
                   {u8"δ", ToMathematica(polynomial2.origin())},
                   {" ", ""},
                   {"\n", ""}}),
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
        absl::StrReplaceAll(
            u8R"(Function[Plus[
                  α,
                  Times[β,Sin[Times[ω,Subtract[#,δ]]]],
                  Times[γ,Cos[Times[ω,Subtract[#,δ]]]]]])",
            {{u8"α", ToMathematicaBody(secular, /*express_in=*/std::nullopt)},
             {u8"β", ToMathematicaBody(sin, /*express_in=*/std::nullopt)},
             {u8"γ", ToMathematicaBody(cos, /*express_in=*/std::nullopt)},
             {u8"δ", ToMathematica(t0)},
             {u8"ω", ToMathematica(4 * Radian / Second)},
             {" ", ""},
             {"\n", ""}}),
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
        absl::StrReplaceAll(
            u8"Function[Piecewise[List[List[α,Between[#,β]]]]]",
            {{u8"α", ToMathematicaBody(series, /*express_in=*/std::nullopt)},
             {u8"β", ToMathematica(std::tuple{interval.min, interval.max})}}),
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
    DiscreteTraject0ry<F> trajectory;
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
