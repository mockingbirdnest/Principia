#include "numerics/polynomial_in_monomial_basis.hpp"

#include <tuple>

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/numerics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/check_well_formedness.hpp"  // 🧙 For PRINCIPIA_CHECK_ILL_FORMED.
#include "testing_utilities/matchers.hpp"

#define PRINCIPIA_USE_IACA 0
#if PRINCIPIA_USE_IACA
#include "intel/iacaMarks.h"
#endif

namespace principia {
namespace numerics {

using ::testing::Eq;
using namespace boost::multiprecision;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_matchers;

class PolynomialInMonomialBasisTest : public ::testing::Test {
 public:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using P2V = PolynomialInMonomialBasis<Displacement<World>, Time, 2>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2>;
  using P2P = PolynomialInMonomialBasis<Position<World>, Instant, 2>;
  using P17 = PolynomialInMonomialBasis<Displacement<World>, Time, 17>;

 protected:
  PolynomialInMonomialBasisTest()
      : coefficients_({
            Displacement<World>({0 * Metre,
                                 0 * Metre,
                                 1 * Metre}),
            Velocity<World>({0 * Metre / Second,
                             1 * Metre / Second,
                             0 * Metre / Second}),
            Vector<Acceleration, World>({1 * Metre / Second / Second,
                                         0 * Metre / Second / Second,
                                         0 * Metre / Second / Second})}) {}

  P2V::Coefficients const coefficients_;
};

#if PRINCIPIA_USE_IACA
// A convenient skeleton for analysing code with IACA.
TEST_F(PolynomialInMonomialBasisTest, DISABLED_IACA) {
  constexpr int degree = 17;
  using E = EstrinEvaluator<Displacement<World>, Time, degree>;
  using P = PolynomialInMonomialBasis<Displacement<World>,
                                      Time,
                                      degree,
                                      EstrinEvaluator>;
  P::Coefficients const coefficients;

  auto iaca = [](P::Coefficients const& c, Time const& t) {
    IACA_VC64_START;
    auto const result = E::Evaluate(c, t);
    IACA_VC64_END;
    return result;
  };
  CHECK_EQ(iaca(coefficients, 2 * Second), iaca(coefficients, 2 * Second));
}
#endif

// Check that coefficients can be accessed and have the right type.
TEST_F(PolynomialInMonomialBasisTest, Coefficients) {
  Displacement<World> const d = std::get<0>(coefficients_);
  Velocity<World> const v = std::get<1>(coefficients_);
  EXPECT_EQ(1 * Metre, d.coordinates().z);
  EXPECT_EQ(1 * Metre / Second, v.coordinates().y);
}

// Check that a polynomial can be constructed and evaluated.
TEST_F(PolynomialInMonomialBasisTest, Evaluate2V) {
  P2V const p(coefficients_, with_evaluator<Horner>);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p(0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial can be for an affine argument.
TEST_F(PolynomialInMonomialBasisTest, Evaluate2A) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2A const p(coefficients_, t0);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p(t0 + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));

  // This compiles.
  p.Primitive();
}

// Check that a polynomial can return an affine value.
TEST_F(PolynomialInMonomialBasisTest, Evaluate2P) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2P const p({World::origin + std::get<0>(coefficients_),
               std::get<1>(coefficients_),
               std::get<2>(coefficients_)},
              t0);
  EXPECT_EQ(2, p.degree());
  Position<World> const d = p(t0 + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(World::origin + Displacement<World>({0.25 * Metre,
                                                                   0.5 * Metre,
                                                                   1 * Metre}),
                              0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

PRINCIPIA_CHECK_WELL_FORMED(
    p.Primitive(),
    with_variable<PolynomialInMonomialBasisTest::P2A> p);
PRINCIPIA_CHECK_ILL_FORMED(p.Primitive(),
                           with_variable<PolynomialInMonomialBasisTest::P2P> p);

// Check that a polynomial of high order may be declared.
TEST_F(PolynomialInMonomialBasisTest, Evaluate17) {
  P17::Coefficients const coefficients;
  P17 const p(coefficients);
  EXPECT_EQ(17, p.degree());
  Displacement<World> const d = p(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}), 0));
}

// Check that a conversion to increase the degree works.
TEST_F(PolynomialInMonomialBasisTest, Conversion) {
  P2V const p2v(coefficients_);
  P17 const p17 = P17(p2v);
  Displacement<World> const d = p17(0.5 * Second);
  Velocity<World> const v = p17.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

TEST_F(PolynomialInMonomialBasisTest, VectorSpace) {
  P2V const p2v(coefficients_);
  {
    auto const p = p2v + p2v;
    auto const actual = p(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v - p2v;
    auto const actual = p(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = 3.0 * Joule * p2v;
    auto const actual = p(0 * Second);
    auto const expected = Vector<Product<Energy, Length>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v * (3.0 * Joule);
    auto const actual = p(0 * Second);
    auto const expected = Vector<Product<Length, Energy>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v / (4.0 * Joule);
    auto const actual = p(0 * Second);
    auto const expected = Vector<Quotient<Length, Energy>, World>(
                  {0 * Metre / Joule, 0 * Metre / Joule, 0.25 * Metre / Joule});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
}

TEST_F(PolynomialInMonomialBasisTest, Ring) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});
  auto const p = p2 * p3;
  {
    auto const actual = p(0 * Second);
    EXPECT_THAT(actual, AlmostEquals(2 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-8 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(-1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-80 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-350 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p(-2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-518 * Ampere * Kelvin, 0));
  }
}

TEST_F(PolynomialInMonomialBasisTest, Affine) {
  using P0A = PolynomialInMonomialBasis<Instant, Time, 0>;
  using P0V = PolynomialInMonomialBasis<Time, Time, 0>;

  P0A const p0a(std::tuple{Instant() + 1 * Second});
  P0V const p0v(std::tuple{2 * Second});
#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
  {
    P0A const p = p0v + Instant();
    EXPECT_THAT(p(3 * Second), AlmostEquals(Instant() + 2 * Second, 0));
  }
#endif
  {
    P0A const p =  Instant() + p0v;
    EXPECT_THAT(p(3 * Second), AlmostEquals(Instant() + 2 * Second, 0));
  }
  {
    P0V const p = p0a - Instant();
    EXPECT_THAT(p(3 * Second), AlmostEquals(1 * Second, 0));
  }
  {
    P0V const p = Instant() - p0a;
    EXPECT_THAT(p(3 * Second), AlmostEquals(-1 * Second, 0));
  }
}

TEST_F(PolynomialInMonomialBasisTest, Monoid) {
  using P1A = PolynomialInMonomialBasis<Current, Temperature, 1>;
  using P2A = PolynomialInMonomialBasis<Temperature, Instant, 2>;
  using P2V = PolynomialInMonomialBasis<Temperature, Time, 2>;
  using P3 = PolynomialInMonomialBasis<Current, Temperature, 3>;
  Instant const t0;
  P1A const p1a({2 * Ampere,
                 -4 * Ampere / Kelvin}, 3 * Kelvin);
  P2A const p2a({1 * Kelvin,
                 3 * Kelvin / Second,
                 -8 * Kelvin / Second / Second}, t0 + 4 * Second);
  P2V const p2v({1 * Kelvin,
                 3 * Kelvin / Second,
                 -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Kelvin,
               3 * Ampere / Kelvin / Kelvin,
               1 * Ampere / Kelvin / Kelvin / Kelvin});
  auto const pa = Compose(p3, p2a);
  auto const pv = Compose(p3, p2v);
  {
    auto const actual_a = pa(t0 + 0 * Second);
    auto const actual_v = pv(0 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-2'627'098 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(2 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 + 1 * Second);
    auto const actual_v = pv(1 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-492'478 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(2 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 - 1 * Second);
    auto const actual_v = pv(-1 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-9'662'098 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-658 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 + 2 * Second);
    auto const actual_v = pv(2 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-46'396 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-13648 * Ampere, 0));
  }
  {
    auto const actual_a = pa(t0 - 2 * Second);
    auto const actual_v = pv(-2 * Second);
    EXPECT_THAT(actual_a, AlmostEquals(-28'092'328 * Ampere, 0));
    EXPECT_THAT(actual_v, AlmostEquals(-46396 * Ampere, 0));
  }
  {
    auto const actual = Compose(p1a, p2a)(t0 + 1 * Second);
    EXPECT_THAT(actual, AlmostEquals(334 * Ampere, 0));
  }
}

TEST_F(PolynomialInMonomialBasisTest, PointwiseInnerProduct) {
  P2V::Coefficients const coefficients({
      Displacement<World>({0 * Metre,
                           2 * Metre,
                           3 * Metre}),
      Velocity<World>({-1 * Metre / Second,
                       1 * Metre / Second,
                       0 * Metre / Second}),
      Vector<Acceleration, World>({1 * Metre / Second / Second,
                                   1 * Metre / Second / Second,
                                   -2 * Metre / Second / Second})});
  P2V const p2va(coefficients_);
  P2V const p2vb(coefficients);

  auto const p = PointwiseInnerProduct(p2va, p2vb);
  {
    auto const actual = p(0 * Second);
    EXPECT_THAT(actual, AlmostEquals(3 * Metre * Metre, 0));
  }
  {
    auto const actual = p(1 * Second);
    EXPECT_THAT(actual, AlmostEquals(5 * Metre * Metre, 0));
  }
  {
    auto const actual = p(-1 * Second);
    EXPECT_THAT(actual, AlmostEquals(1 * Metre * Metre, 0));
  }
  {
    auto const actual = p(2 * Second);
    EXPECT_THAT(actual, AlmostEquals(19 * Metre * Metre, 0));
  }
  {
    auto const actual = p(-2 * Second);
    EXPECT_THAT(actual, AlmostEquals(11 * Metre * Metre, 0));
  }
}

TEST_F(PolynomialInMonomialBasisTest, AtOrigin) {
  Instant const t0 = Instant() + 3 * Second;
  P2A const p(coefficients_, t0);
  P2A const q = p.AtOrigin(Instant() - 2 * Second);
  for (Instant t = Instant() - 10 * Second;
       t < Instant() + 10 * Second;
       t += 0.3 * Second) {
    EXPECT_THAT(q(t), AlmostEquals(p(t), 0, 942));
  }
}

TEST_F(PolynomialInMonomialBasisTest, Derivative) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});

  EXPECT_EQ(3 * Kelvin / Second,
            p2.Derivative<1>()(0 * Second));
  EXPECT_EQ(-16 * Kelvin / Second / Second,
            p2.Derivative<2>()(0 * Second));

  EXPECT_EQ(-4 * Ampere / Second,
            p3.Derivative<1>()(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second,
            p3.Derivative<2>()(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second / Second,
            p3.Derivative<3>()(0 * Second));
}

TEST_F(PolynomialInMonomialBasisTest, PrimitiveIntegrate) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});

  EXPECT_THAT(p2.Primitive()(0 * Second),
              AlmostEquals(0 * Kelvin * Second, 0));
  EXPECT_THAT(p2.Primitive()(1 * Second),
              AlmostEquals(-1.0 / 6.0 * Kelvin * Second, 5));
  EXPECT_THAT(p2.Primitive()(-1 * Second),
              AlmostEquals(19.0 / 6.0 * Kelvin * Second, 1));
  EXPECT_THAT(p2.Primitive()(2 * Second),
              AlmostEquals(-40.0 / 3.0 * Kelvin * Second, 1));

  EXPECT_THAT(p2.Integrate(-1 * Second, 2 * Second),
              AlmostEquals(-99.0 / 6.0 * Kelvin * Second, 3));
}

TEST_F(PolynomialInMonomialBasisTest, EvaluateConstant) {
  PolynomialInMonomialBasis<Entropy, Time, 0> const horner_boltzmann(
      std::make_tuple(BoltzmannConstant), with_evaluator<Horner>);
  PolynomialInMonomialBasis<Entropy, Time, 0> const estrin_boltzmann(
      std::make_tuple(BoltzmannConstant), with_evaluator<Estrin>);
  EXPECT_THAT(horner_boltzmann(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(estrin_boltzmann(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(horner_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
  EXPECT_THAT(estrin_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
}

TEST_F(PolynomialInMonomialBasisTest, EvaluateLinear) {
  PolynomialInMonomialBasis<Length, Time, 1> const
      horner_light({0 * Metre, SpeedOfLight}, with_evaluator<Horner>);
  PolynomialInMonomialBasis<Length, Time, 1> const
      estrin_light({0 * Metre, SpeedOfLight}, with_evaluator<Estrin>);
  constexpr Length light_second = Second * SpeedOfLight;
  EXPECT_THAT(horner_light(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(estrin_light(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(horner_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
  EXPECT_THAT(estrin_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
}

// Check that polynomials may be used with Boost multiprecision types.
TEST_F(PolynomialInMonomialBasisTest, Boost) {
  using P2i = PolynomialInMonomialBasis<cpp_int, cpp_int, 2>;
  P2i const p2i({1, 3, -8});
  EXPECT_EQ(p2i(3), -62);

  using P2r = PolynomialInMonomialBasis<cpp_rational, cpp_rational, 2>;
  P2r const p2r({1, 3, -8});
  EXPECT_EQ(p2r(cpp_rational(1, 3)), cpp_rational(10, 9));

  using P2f = PolynomialInMonomialBasis<cpp_bin_float_50, cpp_bin_float_50, 2>;
  P2f const p2f({1, 3, -8});
  EXPECT_EQ(p2f(cpp_bin_float_50("1.25")), cpp_bin_float_50("-7.75"));
}

// Check that polynomials may be serialized.
TEST_F(PolynomialInMonomialBasisTest, Serialization) {
  {
    P2V p2v(coefficients_);
    serialization::Polynomial message;
    p2v.WriteToMessage(&message);
    EXPECT_EQ(2, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_TRUE(extension.has_quantity());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage(message);
    EXPECT_EQ(2, polynomial_read->degree());
    EXPECT_THAT(
        (*polynomial_read)(0.5 * Second),
        AlmostEquals(
            Displacement<World>({0.25 * Metre, 0.5 * Metre, 1 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
  {
    P2A p2a(coefficients_, Instant());
    serialization::Polynomial message;
    p2a.WriteToMessage(&message);
    EXPECT_EQ(2, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(3, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_TRUE(extension.has_point());
    EXPECT_TRUE(extension.point().has_scalar());

    // Simulate a pre-Καραθεοδωρή compatibility read.
    serialization::Polynomial compatibility_message = message;
    compatibility_message
        .MutableExtension(serialization::PolynomialInMonomialBasis::extension)
        ->clear_evaluator();
    auto const polynomial_read =
        Polynomial<Displacement<World>, Instant>::ReadFromMessage<Horner>(
            compatibility_message);
    EXPECT_EQ(2, polynomial_read->degree());
    *polynomial_read =
        std::move(
            dynamic_cast<
                PolynomialInMonomialBasis<Displacement<World>, Instant, 2>&>(
                *polynomial_read))
            .template WithEvaluator<Horner>();
    EXPECT_THAT(
        (*polynomial_read)(Instant() + 0.5 * Second),
        AlmostEquals(
            Displacement<World>({0.25 * Metre, 0.5 * Metre, 1 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
  {
    P17::Coefficients const coefficients;
    P17 p17(coefficients);
    serialization::Polynomial message;
    p17.WriteToMessage(&message);
    EXPECT_EQ(17, message.degree());
    EXPECT_TRUE(message.HasExtension(
        serialization::PolynomialInMonomialBasis::extension));
    auto const& extension = message.GetExtension(
        serialization::PolynomialInMonomialBasis::extension);
    EXPECT_EQ(18, extension.coefficient_size());
    for (auto const& coefficient : extension.coefficient()) {
      EXPECT_TRUE(coefficient.has_multivector());
    }
    EXPECT_TRUE(extension.has_quantity());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage(message);
    EXPECT_EQ(17, polynomial_read->degree());
    *polynomial_read =
        std::move(
            dynamic_cast<
                PolynomialInMonomialBasis<Displacement<World>, Time, 17>&>(
                *polynomial_read))
            .template WithEvaluator<Estrin>();
    EXPECT_THAT((*polynomial_read)(0.5 * Second),
                AlmostEquals(
                    Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
}

TEST_F(PolynomialInMonomialBasisTest, Output) {
  P2V p2v(coefficients_);
  P2A p2a(coefficients_, Instant());
  P17::Coefficients const coefficients;
  P17 p17(coefficients);
  LOG(ERROR) << p2v;
  LOG(ERROR) << p2a;
  LOG(ERROR) << p17;
}

}  // namespace numerics
}  // namespace principia

#undef PRINCIPIA_USE_IACA
