#include "numerics/fast_fourier_transform.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {
namespace internal_fast_fourier_transform {

using quantities::Sqrt;
using testing_utilities::AlmostEquals;
using ::testing::ElementsAre;

class FastFourierTransformTest : public ::testing::Test {
protected:
  using Complex = std::complex<double>;

  template<typename Container, int size_>
  std::array<std::complex<double>, size_> Coefficients(
      FastFourierTransform<Container, size_> const& fft) {
    return fft.transform_;
  }
};

TEST_F(FastFourierTransformTest, Square) {
  using FFT = FastFourierTransform<std::vector<double>, 8>;
  FFT const transform({1, 1, 1, 1, 0, 0, 0, 0});
  EXPECT_THAT(Coefficients(transform),
              ElementsAre(AlmostEquals(Complex{4}, 0),
                          AlmostEquals(Complex{1, -1 - Sqrt(2)}, 4),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, 1 - Sqrt(2)}, 2),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, Sqrt(2) - 1}, 4),
                          AlmostEquals(Complex{0}, 0),
                          AlmostEquals(Complex{1, 1 + Sqrt(2)}, 3)));
}

TEST_F(FastFourierTransformTest, Sin) {
  using FFT = FastFourierTransform<std::vector<double>, 16>;
  // Sin(x) on [0, 7].
  FFT const transform({+0,
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
                       +0.65698659871878909040});
  EXPECT_THAT(
      Coefficients(transform),
      ElementsAre(
          AlmostEquals(Complex{+0.8462402252352593992601464999}, 0),
          AlmostEquals(Complex{+4.0811719972720017737297705383,
                               -5.7509914354728044020475598497}, 1),
          AlmostEquals(Complex{-1.2157306272208018501429984476,
                               +1.7603441806320655123815748876}, 1),
          AlmostEquals(Complex{-0.7372425382844149484214817488,
                               +0.8380310985140470061957305286}, 9),
          AlmostEquals(Complex{-0.6197133589762012113067702090,
                               +0.5183935643893893471214059230}, 0),
          AlmostEquals(Complex{-0.5726936754458237681916166204,
                               +0.3352695668656918501471058013}, 6),
          AlmostEquals(Complex{-0.5504467338298529673542814016,
                               +0.2045798116680680999699630079}, 6),
          AlmostEquals(Complex{-0.5400094598886346723188351187,
                               +0.0975085343055921838098913009}, 61),
          AlmostEquals(Complex{-0.5369114324878041112477204841}, 2),
          AlmostEquals(Complex{-0.5400094598886346723188351187,
                               -0.0975085343055921838098913009}, 35),
          AlmostEquals(Complex{-0.5504467338298529673542814016,
                               -0.2045798116680680999699630079}, 2),
          AlmostEquals(Complex{-0.5726936754458237681916166204,
                               -0.3352695668656918501471058013}, 6),
          AlmostEquals(Complex{-0.6197133589762012113067702090,
                               -0.5183935643893893471214059230}, 0),
          AlmostEquals(Complex{-0.7372425382844149484214817488,
                               -0.8380310985140470061957305286}, 7),
          AlmostEquals(Complex{-1.2157306272208018501429984476,
                               -1.7603441806320655123815748876}, 2),
          AlmostEquals(Complex{+4.0811719972720017737297705383,
                               +5.7509914354728044020475598497}, 1)));
}

}  // namespace internal_fast_fourier_transform
}  // namespace numerics
}  // namespace principia
