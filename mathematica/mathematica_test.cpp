#include "mathematica/mathematica.hpp"

#include <list>
#include <string>
#include <vector>

#include "absl/strings/str_replace.h"
#include "astronomy/orbital_elements.hpp"
#include "boost/multiprecision/cpp_bin_float.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "geometry/point.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/space.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "gtest/gtest.h"
#include "numerics/double_precision.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/piecewise_poisson_series.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "numerics/polynomial_in_чебышёв_basis.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.

namespace principia {
namespace mathematica {

using namespace boost::multiprecision;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_point;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_space;
using namespace principia::geometry::_symmetric_bilinear_form;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_double_precision;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_piecewise_poisson_series;
using namespace principia::numerics::_poisson_series;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::numerics::_polynomial_in_чебышёв_basis;
using namespace principia::numerics::_unbounded_arrays;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class MathematicaTest : public ::testing::Test {
 protected:
  using F = Frame<struct FTag>;
};

TEST_F(MathematicaTest, ToMathematica) {
  {
    EXPECT_EQ(
        absl::StrReplaceAll(
            "List[α,β,γ]",
            {{"α", ToMathematica(1 * Metre, PreserveUnits)},
             {"β", ToMathematica(2 * Second, PreserveUnits)},
             {"γ", ToMathematica(3 * Metre / Second, PreserveUnits)}}),
        ToMathematica(std::tuple{1 * Metre, 2 * Second, 3 * Metre / Second},
                      PreserveUnits));
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
    EXPECT_EQ("Times[2^^10011010001010110011110001001101010111100110111101111,"
              "Power[2,Subtract[-163,52]]]",
              ToMathematica(0x1.3'4567'89AB'CDEFp-163, std::nullopt, 2));
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
    EXPECT_EQ("Times[16^^1C7005218FCD2A07288057EE16EDA202D8B6BE36DC4,"
              "Power[2,Subtract[-542,168]]]",
              ToMathematica(cpp_bin_float_50("1.23456789e-163")));
    EXPECT_EQ("Times[2^^1110001110000000001010010000110001111110011010010101000"
              "0001110010100010000000010101111110111000010110111011011010001000"
              "00001011011000101101101011111000110110110111000100,"
              "Power[2,Subtract[-542,168]]]",
              ToMathematica(
                  cpp_bin_float_50("1.23456789e-163"), std::nullopt, 2));
    EXPECT_EQ("Times[16^^1000000000000000000000000000000000000000000,"
              "Power[2,Subtract[-1024,168]]]",
              ToMathematica(
                  cpp_bin_float_50(ldexp(cpp_bin_float_50("1"), -1024))));
    EXPECT_EQ("0", ToMathematica(cpp_bin_float_50("0.0")));
    EXPECT_EQ("Minus[0]",
              ToMathematica(
                  -cpp_bin_float_50("0.0")));  // Note that "-0.0" parses as 0.
    EXPECT_EQ("Infinity", ToMathematica(cpp_bin_float_50("inf")));
    EXPECT_EQ("Minus[Infinity]", ToMathematica(cpp_bin_float_50("-inf")));
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
        absl::StrReplaceAll("Quaternion[α,β,γ,δ]",
                            {{"α", ToMathematica(1.0)},
                             {"β", ToMathematica(2.0)},
                             {"γ", ToMathematica(3.0)},
                             {"δ", ToMathematica(-4.0)}}),
        ToMathematica(Quaternion(1.0, R3Element<double>(2.0, 3.0, -4.0))));
  }
  {
    EXPECT_EQ(absl::StrReplaceAll(u8R"(Quantity[α," m s^-1"])",
                                  {{"α", ToMathematica(1.5)}}),
              ToMathematica(1.5 * Metre / Second, PreserveUnits));
  }
  {
    DoublePrecision<double> d(3);
    d += 5e-20;
    EXPECT_EQ(absl::StrReplaceAll(
                  "Plus[α,β]",
                  {{"α", ToMathematica(3.0)},
                   {"β", ToMathematica(5e-20)}}),
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
    EXPECT_EQ(ToMathematica(std::tuple{dof.position(), dof.velocity()},
                            PreserveUnits),
              ToMathematica(dof, PreserveUnits));
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
    DiscreteTrajectory<F> trajectory;
    EXPECT_OK(trajectory.Append(
        Instant(),
        DegreesOfFreedom<F>(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second}))));
    EXPECT_EQ(ToMathematica(std::tuple{trajectory.front().time,
                                       trajectory.front().degrees_of_freedom},
                            PreserveUnits),
              ToMathematica(trajectory.front(), PreserveUnits));
  }
  {
    OrbitalElements::EquinoctialElements elements{
        Instant(), 1 * Metre, 2, 3, 4 * Radian, 5, 6, 7, 8};
    EXPECT_EQ(
        ToMathematica(std::tuple{
            Instant(), 1 * Metre, 2.0, 3.0, 4 * Radian, 5.0, 6.0, 7.0, 8.0},
            PreserveUnits),
        ToMathematica(elements, PreserveUnits));
  }
  {
    PolynomialInMonomialBasis<Length, Time, 2> polynomial1(
        {2 * Metre, -3 * Metre / Second, 4 * Metre / Second / Second});
    EXPECT_EQ(
        absl::StrReplaceAll(
            "Function[Plus[α,Times[β,#],Times[γ,Power[#,2]]]]",
            {{"α", ToMathematica(2 * Metre, PreserveUnits)},
             {"β", ToMathematica(-3 * Metre / Second, PreserveUnits)},
             {"γ", ToMathematica(4 * Metre / Second / Second, PreserveUnits)}}),
        ToMathematica(polynomial1, PreserveUnits));
    PolynomialInMonomialBasis<Length, Instant, 2> polynomial2(
        {5 * Metre, 6 * Metre / Second, -7 * Metre / Second / Second},
        Instant() + 2 * Second);
    EXPECT_EQ(
        absl::StrReplaceAll(
            u8R"(Function[Plus[
                        α,
                        Times[β,Subtract[#,δ]],
                        Times[γ,Power[Subtract[#,δ],2]]]])",
            {{"α", ToMathematica(5 * Metre, PreserveUnits)},
             {"β", ToMathematica(6 * Metre / Second, PreserveUnits)},
             {"γ", ToMathematica(-7 * Metre / Second / Second, PreserveUnits)},
             {"δ", ToMathematica(polynomial2.origin(), PreserveUnits)},
             {" ", ""},
             {"\n", ""}}),
        ToMathematica(polynomial2, PreserveUnits));
  }
  {
    PolynomialInЧебышёвBasis<Length, Time, 2> series(
        {1 * Metre, -3 * Metre, 2 * Metre}, 2 * Second, 4 * Second);
    EXPECT_EQ(absl::StrReplaceAll(
                  R"(Function[Plus[
                            Times[α,ChebyshevT[0,Divide[Subtract[#,δ],ε]]],
                            Times[β,ChebyshevT[1,Divide[Subtract[#,δ],ε]]],
                            Times[γ,ChebyshevT[2,Divide[Subtract[#,δ],ε]]]]])",
                  {{"α", ToMathematica(1 * Metre, PreserveUnits)},
                   {"β", ToMathematica(-3 * Metre, PreserveUnits)},
                   {"γ", ToMathematica(2 * Metre, PreserveUnits)},
                   {"δ", ToMathematica(3 * Second, PreserveUnits)},
                   {"ε", ToMathematica(1 * Second, PreserveUnits)},
                   {" ", ""},
                   {"\n", ""}}),
              ToMathematica(series, PreserveUnits));
  }
  {
    using Series = PoissonSeries<double, 0, 0>;
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
            {{"α", ToMathematicaBody(secular, PreserveUnits)},
             {"β", ToMathematicaBody(sin, PreserveUnits)},
             {"γ", ToMathematicaBody(cos, PreserveUnits)},
             {"δ", ToMathematica(t0, PreserveUnits)},
             {"ω", ToMathematica(4 * Radian / Second, PreserveUnits)},
             {" ", ""},
             {"\n", ""}}),
        ToMathematica(series, PreserveUnits));
  }
  {
    using PiecewiseSeries = PiecewisePoissonSeries<double, 0, 0>;
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
            "Function[Piecewise[List[List[α,Between[#,β]]]]]",
            {{"α", ToMathematicaBody(series, PreserveUnits)},
             {"β",
              ToMathematica(std::tuple{interval.min, interval.max},
                            PreserveUnits)}}),
        ToMathematica(pw, PreserveUnits));
  }
  {
    std::optional<std::string> opt1;
    std::optional<std::string> opt2("foo");
    EXPECT_EQ("List[]", ToMathematica(opt1));
    EXPECT_EQ("List[\"foo\"]", ToMathematica(opt2));
  }
}

TEST_F(MathematicaTest, Rule) {
  EXPECT_EQ("Rule[option," + ToMathematica(3.0) + "]", Rule("option", 3.0));
}

TEST_F(MathematicaTest, Set) {
  EXPECT_EQ("Set[var," + ToMathematica(3.0) + "];\n", Set("var", 3.0));
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
    EXPECT_EQ(
        ToMathematica(3.0),
        ToMathematica(3.0 * Metre / Second / Second, ExpressInSIUnits));
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
    EXPECT_OK(trajectory.Append(
        Instant(),
        DegreesOfFreedom<F>(
            F::origin +
                Displacement<F>({2.0 * Metre, 3.0 * Metre, -4.0 * Metre}),
            Velocity<F>({-1.0 * Metre / Second,
                         -5.0 * Metre / Second,
                         8.0 * Metre / Second}))));
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

}  // namespace mathematica
}  // namespace principia
