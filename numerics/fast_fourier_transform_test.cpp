#include "numerics/fast_fourier_transform.hpp"

#include <algorithm>
#include <random>
#include <vector>

#include "geometry/complexification.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Sqrt;
using quantities::Time;
using quantities::Voltage;
using quantities::si::Metre;
using quantities::si::Second;
using quantities::si::Volt;
using testing_utilities::AlmostEquals;
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::Lt;
using ::testing::Pair;

class FastFourierTransformTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

  using Complex = Complexification<double>;

  template<typename Scalar, std::size_t size_>
  std::array<Complexification<double>, size_> Coefficients(
      FastFourierTransform<Scalar, Instant, size_> const& fft) {
    return fft.transform_;
  }
};

TEST_F(FastFourierTransformTest, Square) {
  FastFourierTransform<double, Instant, 8> const transform(
      {1, 1, 1, 1, 0, 0, 0, 0}, 1 * Second);
  EXPECT_THAT(Coefficients(transform),
              ElementsAre(AlmostEquals(Complex{4}, 0),
                          AlmostEquals(Complex{1, -1 - Sqrt(2)}, 0, 1),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, 1 - Sqrt(2)}, 4, 8),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, Sqrt(2) - 1}, 4),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, 1 + Sqrt(2)}, 0, 2)));
}

TEST_F(FastFourierTransformTest, Sin) {
  // Sin(x) on [0, 7].
  FastFourierTransform<double, Instant, 16> const transform(
      {+0,
       +0.44991188055599964373,
       +0.80360826369441117592,
       +0.98544972998846018066,
       +0.95654873748436662401,
       +0.72308588173832461680,
       +0.33498815015590491954,
       -0.12474816864589884767,
       -0.55780658091328209620,
       -0.87157577241358806002,
       -0.99895491709792831520,
       -0.91270346343588987220,
       -0.63126663787232131146,
       -0.21483085764466499644,
       +0.24754738092257664739,
       +0.65698659871878909040},
      1 * Second);
  EXPECT_THAT(
      Coefficients(transform),
      ElementsAre(
          AlmostEquals(Complex{+0.8462402252352593992601464999}, 0),
          AlmostEquals(Complex{+4.0811719972720017737297705383,
                               -5.7509914354728044020475598497}, 1),
          AlmostEquals(Complex{-1.2157306272208018501429984476,
                               +1.7603441806320655123815748876}, 1),
          AlmostEquals(Complex{-0.7372425382844149484214817488,
                               +0.8380310985140470061957305286}, 6),
          AlmostEquals(Complex{-0.6197133589762012113067702090,
                               +0.5183935643893893471214059230}, 0),
          AlmostEquals(Complex{-0.5726936754458237681916166204,
                               +0.3352695668656918501471058013}, 2),
          AlmostEquals(Complex{-0.5504467338298529673542814016,
                               +0.2045798116680680999699630079}, 6),
          AlmostEquals(Complex{-0.5400094598886346723188351187,
                               +0.0975085343055921838098913009}, 29),
          AlmostEquals(Complex{-0.5369114324878041112477204841}, 2),
          AlmostEquals(Complex{-0.5400094598886346723188351187,
                               -0.0975085343055921838098913009}, 3),
          AlmostEquals(Complex{-0.5504467338298529673542814016,
                               -0.2045798116680680999699630079}, 2, 6),
          AlmostEquals(Complex{-0.5726936754458237681916166204,
                               -0.3352695668656918501471058013}, 6, 8),
          AlmostEquals(Complex{-0.6197133589762012113067702090,
                               -0.5183935643893893471214059230}, 0),
          AlmostEquals(Complex{-0.7372425382844149484214817488,
                               -0.8380310985140470061957305286}, 1, 2),
          AlmostEquals(Complex{-1.2157306272208018501429984476,
                               -1.7603441806320655123815748876}, 1, 2),
          AlmostEquals(Complex{+4.0811719972720017737297705383,
                               +5.7509914354728044020475598497}, 1)));

  EXPECT_THAT(
      transform.PowerSpectrum(),
      ElementsAre(Pair(AlmostEquals(0 * Radian / Second, 0),
                       AlmostEquals(0.7161225188062225589818894003, 0)),
                  Pair(AlmostEquals(π / 8 * Radian / Second, 0),
                       AlmostEquals(49.729867362198687411669694816, 0)),
                  Pair(AlmostEquals(π / 4 * Radian / Second, 0),
                       AlmostEquals(4.5768125922478623650817219388, 1)),
                  Pair(AlmostEquals(3 * π / 8 * Radian / Second, 0),
                       AlmostEquals(1.2458226823327073992355624995, 3, 4)),
                  Pair(AlmostEquals(π / 2 * Radian / Second, 0),
                       AlmostEquals(0.6527765348939019854655626516, 1)),
                  Pair(AlmostEquals(5 * π / 8 * Radian / Second, 0),
                       AlmostEquals(0.4403837283619551481413056610, 2, 3)),
                  Pair(AlmostEquals(3 * π / 4 * Radian / Second, 0),
                       AlmostEquals(0.3448445061260952118899789115, 1, 5)),
                  Pair(AlmostEquals(7 * π / 8 * Radian / Second, 0),
                       AlmostEquals(0.3011181310316397868684530905, 10, 18)),
                  Pair(AlmostEquals(π * Radian / Second, 0),
                       AlmostEquals(0.2882738863361058320489546747, 4)),
                  Pair(AlmostEquals(9 * π / 8 * Radian / Second, 0),
                       AlmostEquals(0.3011181310316397868684530905, 0)),
                  Pair(AlmostEquals(5 * π / 4 * Radian / Second, 0),
                       AlmostEquals(0.3448445061260952118899789115, 1, 4)),
                  Pair(AlmostEquals(11 * π / 8 * Radian / Second, 0),
                       AlmostEquals(0.4403837283619551481413056610, 4, 15)),
                  Pair(AlmostEquals(3 * π / 2 * Radian / Second, 0),
                       AlmostEquals(0.6527765348939019854655626516, 1)),
                  Pair(AlmostEquals(13 * π / 8 * Radian / Second, 0),
                       AlmostEquals(1.2458226823327073992355624995, 1)),
                  Pair(AlmostEquals(7 * π / 4 * Radian / Second, 0),
                       AlmostEquals(4.5768125922478623650817219388, 1, 2)),
                  Pair(AlmostEquals(15 * π / 8 * Radian / Second, 0),
                       AlmostEquals(49.729867362198687411669694816, 0, 2))));
}

TEST_F(FastFourierTransformTest, Mode) {
  using FFT = FastFourierTransform<Length, Instant, 1 << 16>;
  AngularFrequency const ω = 666 * π / FFT::size * Radian / Second;
  Time const Δt = 1 * Second;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> noise(-0.5, 0.5);
  std::vector<Length> signal;
  for (int n = 0; n < FFT::size; ++n) {
    signal.push_back((Sin(n * ω * Δt) + noise(random)) * Metre);
  }

  // Won't fit on the stack.
  auto transform = std::make_unique<FFT>(signal, Δt);

  auto const mode = transform->Mode();
  EXPECT_THAT(mode.midpoint(), AlmostEquals(ω, 0));
  EXPECT_THAT(mode.measure(),
              AlmostEquals(4 * π / FFT::size * Radian / Second, 24));
}

TEST_F(FastFourierTransformTest, Vector) {
  using FFT = FastFourierTransform<Displacement<World>, Instant, 1 << 16>;
  AngularFrequency const ω = 666 * π / FFT::size * Radian / Second;
  Time const Δt = 1 * Second;
  std::vector<Displacement<World>> signal;
  for (int n = 0; n < FFT::size; ++n) {
    signal.push_back(Displacement<World>({Sin(n * ω * Δt) * Metre,
                                          Cos(n * ω * Δt) * Metre,
                                          Sin(2 * n * ω * Δt) * Metre}));
  }

  // Won't fit on the stack.
  auto transform = std::make_unique<FFT>(signal, Δt);

  auto const mode = transform->Mode();
  EXPECT_THAT(mode.midpoint(), AlmostEquals(ω, 0));
  EXPECT_THAT(mode.measure(),
              AlmostEquals(4 * π / FFT::size * Radian / Second, 24));
}

TEST_F(FastFourierTransformTest, Inverse) {
  constexpr int n = 4;
  constexpr Time Δt = 1 * Second;
  std::array<Voltage, n> v{{1 * Volt, 0 * Volt, 1 * Volt, 0 * Volt}};
  FastFourierTransform<Voltage, Instant, n> const V(v, Δt);

  FastFourierTransform<Voltage, AngularFrequency, n> const nv(
      {V[0].real_part(), V[1].real_part(), V[2].real_part(), V[3].real_part()},
      V.frequency(1) - V.frequency(0));
  EXPECT_THAT((std::array{nv[0].real_part() / n,
                          nv[1].real_part() / n,
                          nv[2].real_part() / n,
                          nv[3].real_part() / n}),
              ElementsAreArray(v));
  EXPECT_THAT(nv.frequency(1) - nv.frequency(0), AlmostEquals(Δt, 0));
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
