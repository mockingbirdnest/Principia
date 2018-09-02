
#include "numerics/combinatorics.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

TEST(Combinatorics, Binomial) {
  // Pascal's triangle.
  static_assert(Binomial(0, 0) == 1, "");

  static_assert(Binomial(1, 0) == 1, "");
  static_assert(Binomial(1, 1) == 1, "");

  static_assert(Binomial(2, 0) == 1, "");
  static_assert(Binomial(2, 1) == 2, "");
  static_assert(Binomial(2, 2) == 1, "");

  static_assert(Binomial(3, 0) == 1, "");
  static_assert(Binomial(3, 1) == 3, "");
  static_assert(Binomial(3, 2) == 3, "");
  static_assert(Binomial(3, 3) == 1, "");

  static_assert(Binomial(4, 0) == 1, "");
  static_assert(Binomial(4, 1) == 4, "");
  static_assert(Binomial(4, 2) == 6, "");
  static_assert(Binomial(4, 3) == 4, "");
  static_assert(Binomial(4, 4) == 1, "");
}

TEST(Combinatorics, DoubleFactorial) {
  // https://oeis.org/A000165, https://oeis.org/A001147.
  static_assert(DoubleFactorial(0) == 1, "");
  static_assert(DoubleFactorial(1) == 1, "");
  static_assert(DoubleFactorial(2) == 2, "");
  static_assert(DoubleFactorial(3) == 3, "");
  static_assert(DoubleFactorial(4) == 8, "");
  static_assert(DoubleFactorial(5) == 15, "");
  static_assert(DoubleFactorial(6) == 48, "");
  static_assert(DoubleFactorial(7) == 105, "");
  static_assert(DoubleFactorial(8) == 384, "");
  static_assert(DoubleFactorial(9) == 945, "");
  static_assert(DoubleFactorial(10) == 3'840, "");
  static_assert(DoubleFactorial(11) == 10'395, "");
  static_assert(DoubleFactorial(12) == 46'080, "");
  static_assert(DoubleFactorial(13) == 135'135, "");
  static_assert(DoubleFactorial(14) == 645'120, "");
  static_assert(DoubleFactorial(15) == 2'027'025, "");
  static_assert(DoubleFactorial(16) == 10'321'920, "");
  static_assert(DoubleFactorial(17) == 34'459'425, "");
  static_assert(DoubleFactorial(18) == 185'794'560, "");
  static_assert(DoubleFactorial(19) == 654'729'075, "");
  static_assert(DoubleFactorial(20) == 3'715'891'200, "");
  static_assert(DoubleFactorial(21) == 13'749'310'575, "");
  static_assert(DoubleFactorial(22) == 81'749'606'400, "");
  static_assert(DoubleFactorial(23) == 316'234'143'225, "");
  static_assert(DoubleFactorial(24) == 1'961'990'553'600, "");
  static_assert(DoubleFactorial(25) == 7'905'853'580'625, "");
  static_assert(DoubleFactorial(26) == 51'011'754'393'600, "");
  static_assert(DoubleFactorial(27) == 213'458'046'676'875, "");
  static_assert(DoubleFactorial(28) == 1'428'329'123'020'800, "");
  static_assert(DoubleFactorial(29) == 6'190'283'353'629'375, "");
  static_assert(DoubleFactorial(30) == 42'849'873'690'624'000, "");
  static_assert(DoubleFactorial(31) == 191'898'783'962'510'625, "");
  static_assert(DoubleFactorial(32) == 1'371'195'958'099'968'000, "");
  static_assert(DoubleFactorial(33) == 6'332'659'870'762'850'625, "");
}

TEST(Combinatorics, Factorial) {
  // https://oeis.org/A000142.
  static_assert(Factorial(0) == 1, "");
  static_assert(Factorial(1) == 1, "");
  static_assert(Factorial(2) == 2, "");
  static_assert(Factorial(3) == 6, "");
  static_assert(Factorial(4) == 24, "");
  static_assert(Factorial(5) == 120, "");
  static_assert(Factorial(6) == 720, "");
  static_assert(Factorial(7) == 5'040, "");
  static_assert(Factorial(8) == 40'320, "");
  static_assert(Factorial(9) == 362'880, "");
  static_assert(Factorial(10) == 3'628'800, "");
  static_assert(Factorial(11) == 39'916'800, "");
  static_assert(Factorial(12) == 479'001'600, "");
  static_assert(Factorial(13) == 6'227'020'800, "");
  static_assert(Factorial(14) == 87'178'291'200, "");
  static_assert(Factorial(15) == 1'307'674'368'000, "");
  static_assert(Factorial(16) == 20'922'789'888'000, "");
  static_assert(Factorial(17) == 355'687'428'096'000, "");
  static_assert(Factorial(18) == 6'402'373'705'728'000, "");
  static_assert(Factorial(19) == 121'645'100'408'832'000, "");
  static_assert(Factorial(20) == 2'432'902'008'176'640'000, "");
}

TEST(Combinatorics, FallingFactorial) {
  // http://oeis.org/A068424.
  static_assert(FallingFactorial(0, 0) == 1, "");
  static_assert(FallingFactorial(1, 0) == 1, "");
  static_assert(FallingFactorial(1, 1) == 1, "");
  static_assert(FallingFactorial(2, 0) == 1, "");
  static_assert(FallingFactorial(2, 1) == 2, "");
  static_assert(FallingFactorial(2, 2) == 2, "");
  static_assert(FallingFactorial(3, 0) == 1, "");
  static_assert(FallingFactorial(3, 1) == 3, "");
  static_assert(FallingFactorial(3, 2) == 6, "");
  static_assert(FallingFactorial(3, 3) == 6, "");
  static_assert(FallingFactorial(4, 0) == 1, "");
  static_assert(FallingFactorial(4, 1) == 4, "");
  static_assert(FallingFactorial(4, 2) == 12, "");
  static_assert(FallingFactorial(4, 3) == 24, "");
  static_assert(FallingFactorial(4, 4) == 24, "");
  static_assert(FallingFactorial(5, 0) == 1, "");
  static_assert(FallingFactorial(5, 1) == 5, "");
  static_assert(FallingFactorial(5, 2) == 20, "");
  static_assert(FallingFactorial(5, 3) == 60, "");
  static_assert(FallingFactorial(5, 4) == 120, "");
  static_assert(FallingFactorial(5, 5) == 120, "");
}

}  // namespace numerics
}  // namespace principia
