
#include "numerics/polynomial.hpp"

#include <tuple>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/numerics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"

#define PRINCIPIA_USE_IACA 0
#if PRINCIPIA_USE_IACA
#include "intel/iacaMarks.h"
#endif

namespace principia {

using geometry::Frame;
using geometry::Displacement;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Vector;
using geometry::Velocity;
using quantities::Acceleration;
using quantities::Energy;
using quantities::Entropy;
using quantities::Length;
using quantities::Product;
using quantities::Quotient;
using quantities::Current;
using quantities::Temperature;
using quantities::Time;
using quantities::constants::BoltzmannConstant;
using quantities::constants::SpeedOfLight;
using quantities::si::Ampere;
using quantities::si::Joule;
using quantities::si::Kelvin;
using quantities::si::Metre;
using quantities::si::Second;
using quantities::si::Watt;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;
using ::testing::Eq;

namespace numerics {

class PolynomialTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using P2V = PolynomialInMonomialBasis<Displacement<World>, Time, 2,
                                        HornerEvaluator>;
  using P2A = PolynomialInMonomialBasis<Displacement<World>, Instant, 2,
                                        HornerEvaluator>;
  using P17 = PolynomialInMonomialBasis<Displacement<World>, Time, 17,
                                        EstrinEvaluator>;

  PolynomialTest()
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
TEST_F(PolynomialTest, DISABLED_IACA) {
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
TEST_F(PolynomialTest, Coefficients) {
  Displacement<World> const d = std::get<0>(coefficients_);
  Velocity<World> const v = std::get<1>(coefficients_);
  EXPECT_EQ(1 * Metre, d.coordinates().z);
  EXPECT_EQ(1 * Metre / Second, v.coordinates().y);
}

// Check that a polynomial can be constructed and evaluated.
TEST_F(PolynomialTest, Evaluate2V) {
  P2V const p(coefficients_);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p.Evaluate(0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial can be for an affine argument.
TEST_F(PolynomialTest, Evaluate2A) {
  Instant const t0 = Instant() + 0.3 * Second;
  P2A const p(coefficients_, t0);
  EXPECT_EQ(2, p.degree());
  Displacement<World> const d = p.Evaluate(t0 + 0.5 * Second);
  Velocity<World> const v = p.EvaluateDerivative(t0 + 0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0.25 * Metre,
                                                   0.5 * Metre,
                                                   1 * Metre}), 0));
  EXPECT_THAT(v, AlmostEquals(Velocity<World>({1 * Metre / Second,
                                               1 * Metre / Second,
                                               0 * Metre / Second}), 0));
}

// Check that a polynomial of high order may be declared.
TEST_F(PolynomialTest, Evaluate17) {
  P17::Coefficients const coefficients;
  P17 const p(coefficients);
  EXPECT_EQ(17, p.degree());
  Displacement<World> const d = p.Evaluate(0.5 * Second);
  EXPECT_THAT(d, AlmostEquals(Displacement<World>({0 * Metre,
                                                   0 * Metre,
                                                   0 * Metre}), 0));
}

TEST_F(PolynomialTest, VectorSpace) {
  P2V const p2v(coefficients_);
  {
    auto const p = p2v + p2v;
    auto const actual = p.Evaluate(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v - p2v;
    auto const actual = p.Evaluate(0 * Second);
    auto const expected =
        Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = 3.0 * Joule * p2v;
    auto const actual = p.Evaluate(0 * Second);
    auto const expected = Vector<Product<Energy, Length>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v * (3.0 * Joule);
    auto const actual = p.Evaluate(0 * Second);
    auto const expected = Vector<Product<Length, Energy>, World>(
                  {0 * Joule * Metre, 0 * Joule * Metre, 3 * Joule * Metre});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
  {
    auto const p = p2v / (4.0 * Joule);
    auto const actual = p.Evaluate(0 * Second);
    auto const expected = Vector<Quotient<Length, Energy>, World>(
                  {0 * Metre / Joule, 0 * Metre / Joule, 0.25 * Metre / Joule});
    EXPECT_THAT(actual, AlmostEquals(expected, 0));
  }
}

TEST_F(PolynomialTest, Ring) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3, HornerEvaluator>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});
  auto const p = p2 * p3;
  {
    auto const actual = p.Evaluate(0 * Second);
    EXPECT_THAT(actual, AlmostEquals(2 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p.Evaluate(1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-8 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p.Evaluate(-1 * Second);
    EXPECT_THAT(actual, AlmostEquals(-80 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p.Evaluate(2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-350 * Ampere * Kelvin, 0));
  }
  {
    auto const actual = p.Evaluate(-2 * Second);
    EXPECT_THAT(actual, AlmostEquals(-518 * Ampere * Kelvin, 0));
  }
}

TEST_F(PolynomialTest, Derivative) {
  using P2 = PolynomialInMonomialBasis<Temperature, Time, 2, HornerEvaluator>;
  using P3 = PolynomialInMonomialBasis<Current, Time, 3, HornerEvaluator>;
  P2 const p2({1 * Kelvin, 3 * Kelvin / Second, -8 * Kelvin / Second / Second});
  P3 const p3({2 * Ampere,
               -4 * Ampere / Second,
               3 * Ampere / Second / Second,
               1 * Ampere / Second / Second / Second});

  EXPECT_EQ(3 * Kelvin / Second,
            p2.Derivative<1>().Evaluate(0 * Second));
  EXPECT_EQ(-16 * Kelvin / Second / Second,
            p2.Derivative<2>().Evaluate(0 * Second));

  EXPECT_EQ(-4 * Ampere / Second,
            p3.Derivative<1>().Evaluate(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second,
            p3.Derivative<2>().Evaluate(0 * Second));
  EXPECT_EQ(6 * Ampere / Second / Second / Second,
            p3.Derivative<3>().Evaluate(0 * Second));
}

TEST_F(PolynomialTest, EvaluateConstant) {
  PolynomialInMonomialBasis<Entropy, Time, 0, HornerEvaluator> const
      horner_boltzmann(std::make_tuple(BoltzmannConstant));
  PolynomialInMonomialBasis<Entropy, Time, 0, EstrinEvaluator> const
      estrin_boltzmann(std::make_tuple(BoltzmannConstant));
  EXPECT_THAT(horner_boltzmann.Evaluate(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(estrin_boltzmann.Evaluate(1729 * Second), Eq(BoltzmannConstant));
  EXPECT_THAT(horner_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
  EXPECT_THAT(estrin_boltzmann.EvaluateDerivative(1729 * Second),
              Eq(0 * Watt / Kelvin));
}

TEST_F(PolynomialTest, EvaluateLinear) {
  PolynomialInMonomialBasis<Length, Time, 1, HornerEvaluator> const
      horner_light({0 * Metre, SpeedOfLight});
  PolynomialInMonomialBasis<Length, Time, 1, EstrinEvaluator> const
      estrin_light({0 * Metre, SpeedOfLight});
  constexpr Length light_second = Second * SpeedOfLight;
  EXPECT_THAT(horner_light.Evaluate(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(estrin_light.Evaluate(1729 * Second), Eq(1729 * light_second));
  EXPECT_THAT(horner_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
  EXPECT_THAT(estrin_light.EvaluateDerivative(1729 * Second), Eq(SpeedOfLight));
}

// Check that polynomials may be serialized.
TEST_F(PolynomialTest, Serialization) {
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
    EXPECT_FALSE(extension.has_origin());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage<HornerEvaluator>(
            message);
    EXPECT_EQ(2, polynomial_read->degree());
    EXPECT_THAT(
        polynomial_read->Evaluate(0.5 * Second),
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
    EXPECT_TRUE(extension.has_origin());
    EXPECT_TRUE(extension.origin().has_scalar());

    auto const polynomial_read =
        Polynomial<Displacement<World>,
                   Instant>::ReadFromMessage<HornerEvaluator>(message);
    EXPECT_EQ(2, polynomial_read->degree());
    EXPECT_THAT(
        polynomial_read->Evaluate(Instant() + 0.5 * Second),
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
    EXPECT_FALSE(extension.has_origin());

    auto const polynomial_read =
        Polynomial<Displacement<World>, Time>::ReadFromMessage<HornerEvaluator>(
            message);
    EXPECT_EQ(17, polynomial_read->degree());
    EXPECT_THAT(polynomial_read->Evaluate(0.5 * Second),
                AlmostEquals(
                    Displacement<World>({0 * Metre, 0 * Metre, 0 * Metre}), 0));
    serialization::Polynomial message2;
    polynomial_read->WriteToMessage(&message2);
    EXPECT_THAT(message2, EqualsProto(message));
  }
}

}  // namespace numerics
}  // namespace principia

#undef PRINCIPIA_USE_IACA
