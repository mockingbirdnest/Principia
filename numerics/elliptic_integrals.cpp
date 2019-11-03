
#include "glog/logging.h"

#include <limits>
#include <tuple>
#include <utility>

#include "numerics/combinatorics.hpp"
#include "numerics/elliptic_integrals.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

// The implementation in this file is derived from [Fuk18] (license: MIT). The
// original code has been translated into C++ and adapted to the needs of this
// project.

namespace principia {
namespace numerics {
namespace internal_elliptic_integrals {

using base::uninitialized;
using quantities::Abs;
using quantities::ArcSin;
using quantities::Angle;
using quantities::ArcTan;
using quantities::ArcTanh;
using quantities::Cos;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Tan;
using quantities::si::Radian;

namespace {

// The functions that compute elliptic integrals of multiple kinds may be called
// without computing the integral of the third kind.  Passing |unused| bypasses
// the corresponding code path at compile time, passing a reference to an
// |Angle| performs the computation.

struct UnusedResult {
  constexpr UnusedResult(base::uninitialized_t) {}
};

inline constexpr UnusedResult unused{uninitialized};

template<typename T>
using EnableIfAngleResult =
    std::enable_if_t<std::is_same_v<T, UnusedResult const> ||
                     std::is_same_v<T, Angle>>;

template<typename T, typename = EnableIfAngleResult<T>>
inline constexpr bool should_compute = !std::is_same_v<T, UnusedResult const>;

// Bulirsch's cel function, [Bul69], [OLBC10], 19.2(iii).
Angle BulirschCel(double kc, double nc, double a, double b);

// Jacobi's nome approximated by a series of the given degree.
template<int degree>
double EllipticNomeQ(double mc);

// Computes Emde’s complete elliptic integrals of the second kind B(m) and D(m),
// where m = 1 - mc.  The method is similar to that described in [Fuk11a].
void FukushimaEllipticBD(double mc, Angle& B_m, Angle& D_m);

// Computes Emde’s complete elliptic integrals of the second kind B(m) and D(m),
// as well as Fukushima’s complete elliptic integral of the third kind J(n, m),
// where m = 1 - mc and n = 1 - nc.  The methods are similar to those described
// in [Fuk11a] and [Fuk12].
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void FukushimaEllipticBDJ(double nc,
                          double mc,
                          Angle& B_m,
                          Angle& D_m,
                          ThirdKind& J_n_m);

// Computes Fukushima's incomplete integrals of the second kind and third kind
// from the cosine of the amplitude: Bc(c|m) = B(arccos c|m),
// Dc(c|m) = D(arccos c|m), Jc(c, n|m) = J(arccos c, n|m), where m = 1 - mc.
// These functions are defined in [Fuk11b], equations (9) and (10), and
// [Fuk12], equation (24).
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void FukushimaEllipticBcDcJc(double c,
                             double n,
                             double mc,
                             Angle& Bc_cǀm,
                             Angle& Dc_cǀm,
                             ThirdKind& Jc_c_nǀm);

// Computes Fukushima's incomplete integrals of the second kind and third kind
// from the sine of the amplitude: Bs(s|m) = B(arcsin s|m),
// Ds(s|m) = D(arcsin s|m), Js(s, n|m) = J(arcsin s, n|m), where m = 1 - mc.
// These functions are defined in [Fuk11b], equations (9) and (10), and
// [Fuk12], equation (24).
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void FukushimaEllipticBsDsJs(double s,
                             double n,
                             double mc,
                             Angle& Bs_sǀm,
                             Angle& Ds_sǀm,
                             ThirdKind& Js_s_nǀm);

// Computes the Maclaurin series expansion ∑ Bₗ(m) yˡ and ∑ Dₗ(m) yˡ used in the
// computation of of Bs and Ds, see [Fuk11b], equation (15).
void FukushimaEllipticBsDsMaclaurinSeries(double y,
                                          double m,
                                          Angle& Σ_Bₗ_m_yˡ,
                                          Angle& Σ_Dₗ_m_yˡ);

// Maclaurin series expansion of Js [Fuk12].
Angle FukushimaEllipticJsMaclaurinSeries(double y, double n, double m);

// Fukushima's T function [Fuk12].
Angle FukushimaT(double t, double h);

// Argument reduction: angle = fractional_part + integer_part * π where
// fractional_part is in [-π/2, π/2].
void Reduce(Angle const& angle,
            Angle& fractional_part,
            std::int64_t& integer_part);

// The common implementation underlying the functions |FukushimaEllipticBDJ| and
// |FukushimaEllipticBD| declared in the header file.
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void FukushimaEllipticBDJ(Angle const& φ,
                          double const n,
                          double const mc,
                          Angle& B_φǀm,
                          Angle& D_φǀm,
                          ThirdKind& J_φ_nǀm);

// Implementation of the B, D, J functions with all arguments reduced.
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void FukushimaEllipticBDJReduced(Angle const& φ,
                                 double const n,
                                 double const mc,
                                 Angle& B_φǀm,
                                 Angle& D_φǀm,
                                 ThirdKind& J_φ_nǀm);

// The common implementation underlying the functions |EllipticFEΠ| and
// |EllipticFE| declared in the header file.
template<typename ThirdKind, typename = EnableIfAngleResult<ThirdKind>>
void EllipticFEΠ(Angle const& φ,
                 double const n,
                 double const mc,
                 Angle& F_φǀm,
                 Angle& E_φǀm,
                 ThirdKind& Π_φ_nǀm);

// A generator for the Maclaurin series for q(m) / m where q is Jacobi's nome
// function.
template<int n, template<typename, typename, int> class Evaluator>
class EllipticNomeQMaclaurin {
  static constexpr auto full_series = std::make_tuple(
    1.0 / 16.0,
    1.0 / 32.0,
    21.0 / 1024.0,
    31.0 / 2048.0,
    6257.0 / 524288.0,
    10293.0 / 1048576.0,
    279025.0 / 33554432.0,
    483127.0 / 67108864.0,
    435506703.0 / 68719476736.0,
    776957575.0 / 137438953472.0,
    22417045555.0 / 4398046511104.0,
    40784671953.0 / 8796093022208.0,
    9569130097211.0 / 2251799813685248.0,
    17652604545791.0 / 4503599627370496.0,
    523910972020563.0 / 144115188075855872.0,
    976501268709949.0 / 288230376151711744.0);

  template<typename>
  struct Generator;

  template<int... k>
  struct Generator<std::index_sequence<k...>> {
    static auto constexpr series = std::make_tuple(std::get<k>(full_series)...);
  };

 public:
  static inline PolynomialInMonomialBasis<double, double, n, Evaluator> const
      polynomial{Generator<std::make_index_sequence<n + 1>>::series};
};

// A generator for the Maclaurin series for Fukushima's elliptic Fs integral.
template<int n, template<typename, typename, int> class Evaluator>
class FukushimaEllipticFsMaclaurin {
  template<typename>
  struct Generator;

  template<int... k>
  struct Generator<std::index_sequence<k...>> {
    static auto constexpr series = std::make_tuple(
        static_cast<double>(DoubleFactorial(2 * k - 1) *
                            DoubleFactorial(2 * n - 2 * k - 1)) /
        static_cast<double>((1 << n) * (2 * n + 1) *
                            Factorial(k) * Factorial(n - k))...);
  };

 public:
  static inline PolynomialInMonomialBasis<double, double, n, Evaluator> const
      polynomial{Generator<std::make_index_sequence<n + 1>>::series};
};

using FukushimaEllipticFsMaclaurin1 =
    FukushimaEllipticFsMaclaurin<1, HornerEvaluator>;
using FukushimaEllipticFsMaclaurin2 =
    FukushimaEllipticFsMaclaurin<2, HornerEvaluator>;
using FukushimaEllipticFsMaclaurin3 =
    FukushimaEllipticFsMaclaurin<3, HornerEvaluator>;
using FukushimaEllipticFsMaclaurin4 =
    FukushimaEllipticFsMaclaurin<4, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin5 =
    FukushimaEllipticFsMaclaurin<5, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin6 =
    FukushimaEllipticFsMaclaurin<6, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin7 =
    FukushimaEllipticFsMaclaurin<7, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin8 =
    FukushimaEllipticFsMaclaurin<8, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin9 =
    FukushimaEllipticFsMaclaurin<9, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin10 =
    FukushimaEllipticFsMaclaurin<10, EstrinEvaluator>;
using FukushimaEllipticFsMaclaurin11 =
    FukushimaEllipticFsMaclaurin<11, EstrinEvaluator>;

// A helper class for building the Maclaurin polynomials for the multivariate Bs
// and Ds.
template<template<typename, typename, int> class Evaluator>
class FukushimaEllipticDsBsMaclaurin {
  template<typename, typename>
  struct Generator;

  template<typename Tuple, int... k>
  struct Generator<Tuple, std::index_sequence<k...>> {
    template<int j>
    static double ComputeBsCoefficient(Tuple const& tuple) {
      if constexpr (j == 0) {
        return 1.0;
      } else {
        return std::get<j>(tuple) -
               std::get<j - 1>(tuple) * (static_cast<double>(2 * j - 1) /
                                         static_cast<double>(2 * j + 1));
      }
    }

    static Tuple ComputeBsCoefficients(Tuple const& tuple) {
      return std::make_tuple(ComputeBsCoefficient<k>(tuple)...);
    }

    static Tuple ComputeDsCoefficients(Tuple const& tuple) {
      return std::make_tuple(
          std::get<k>(tuple) *
          (static_cast<double>(2 * k + 1) / static_cast<double>(2 * k + 3))...);
    }
  };

 public:
  template<typename... Args, int n = sizeof...(Args)>
  static PolynomialInMonomialBasis<double, double, n - 1, Evaluator> const
  MakeBsPolynomial(Args... args) {
    using Tuple = std::tuple<Args...>;
    return PolynomialInMonomialBasis<double, double, n - 1, Evaluator>(
        Generator<Tuple, std::make_index_sequence<n>>::ComputeBsCoefficients(
            std::make_tuple(args...)));
  }

  template<typename... Args, int n = sizeof...(Args)>
  static PolynomialInMonomialBasis<double, double, n - 1, Evaluator> const
  MakeDsPolynomial(Args... args) {
    using Tuple = std::tuple<Args...>;
    return PolynomialInMonomialBasis<double, double, n - 1, Evaluator>(
        Generator<Tuple, std::make_index_sequence<n>>::ComputeDsCoefficients(
            std::make_tuple(args...)));
  }
};

// Maclaurin series for Fukushima's elliptic integral Js.  These are polynomials
// of n that will be used as coefficients of a polynomial in m.  The first index
// is the total degree (l - 1 in Fukushima's notation), the second the degree in
// m (k in Fukushima's notation).
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_0_0(std::make_tuple(1.0 / 3.0));

PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_1_0(std::make_tuple(1.0 / 10.0,
                                                          2.0 / 10.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_1_1(std::make_tuple(1.0 / 10.0));

PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_2_0(std::make_tuple(3.0 / 56.0,
                                                          4.0 / 56.0,
                                                          8.0 / 56.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_2_1(std::make_tuple(2.0 / 56.0,
                                                          4.0 / 56.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_2_2(std::make_tuple(3.0 / 56.0));

PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_3_0(std::make_tuple(5.0 / 144.0,
                                                          6.0 / 144.0,
                                                          8.0 / 144.0,
                                                          16.0 / 144.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_3_1(std::make_tuple(3.0 / 144.0,
                                                          4.0 / 144.0,
                                                          8.0 / 144.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_3_2(std::make_tuple(3.0 / 144.0,
                                                          6.0 / 144.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_3_3(std::make_tuple(5.0 / 144.0));

PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_4_0(std::make_tuple(35.0 / 1408.0,
                                                          40.0 / 1408.0,
                                                          48.0 / 1408.0,
                                                          64.0 / 1408.0,
                                                          128.0 / 1408.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_4_1(std::make_tuple(20.0 / 1408.0,
                                                          24.0 / 1408.0,
                                                          32.0 / 1408.0,
                                                          64.0 / 1408.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_4_2(std::make_tuple(18.0 / 1408.0,
                                                          24.0 / 1408.0,
                                                          48.0 / 1408.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_4_3(std::make_tuple(20.0 / 1408.0,
                                                          40.0 / 1408.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_4_4(std::make_tuple(35.0 / 1408.0));

PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_0(std::make_tuple(63.0 / 3328.0,
                                                          70.0 / 3328.0,
                                                          80.0 / 3328.0,
                                                          96.0 / 3328.0,
                                                          128.0 / 3328.0,
                                                          256.0 / 3328.0));
PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_1(std::make_tuple(35.0 / 3328.0,
                                                          40.0 / 3328.0,
                                                          48.0 / 3328.0,
                                                          64.0 / 3328.0,
                                                          128.0 / 3328.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_2(std::make_tuple(30.0 / 3328.0,
                                                          36.0 / 3328.0,
                                                          48.0 / 3328.0,
                                                          96.0 / 3328.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_3(std::make_tuple(30.0 / 3328.0,
                                                          40.0 / 3328.0,
                                                          80.0 / 3328.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_4(std::make_tuple(35.0 / 3328.0,
                                                          70.0 / 3328.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_5_5(std::make_tuple(63.0 / 3328.0));

PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_0(std::make_tuple(231.0 / 15360.0,
                                                          252.0 / 15360.0,
                                                          280.0 / 15360.0,
                                                          320.0 / 15360.0,
                                                          384.0 / 15360.0,
                                                          512.0 / 15360.0,
                                                          1024.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_1(std::make_tuple(126.0 / 15360.0,
                                                          140.0 / 15360.0,
                                                          160.0 / 15360.0,
                                                          192.0 / 15360.0,
                                                          256.0 / 15360.0,
                                                          512.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_2(std::make_tuple(105.0 / 15360.0,
                                                          120.0 / 15360.0,
                                                          144.0 / 15360.0,
                                                          192.0 / 15360.0,
                                                          384.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_3(std::make_tuple(100.0 / 15360.0,
                                                          120.0 / 15360.0,
                                                          160.0 / 15360.0,
                                                          320.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_4(std::make_tuple(105.0 / 15360.0,
                                                          140.0 / 15360.0,
                                                          280.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_5(std::make_tuple(126.0 / 15360.0,
                                                          252.0 / 15360.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_6_6(std::make_tuple(231.0 / 15360.0));

PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_0(std::make_tuple(429.0 / 34816.0,
                                                          462.0 / 34816.0,
                                                          504.0 / 34816.0,
                                                          560.0 / 34816.0,
                                                          640.0 / 34816.0,
                                                          768.0 / 34816.0,
                                                          1024.0 / 34816.0,
                                                          2048.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_1(std::make_tuple(231.0 / 34816.0,
                                                          252.0 / 34816.0,
                                                          280.0 / 34816.0,
                                                          320.0 / 34816.0,
                                                          384.0 / 34816.0,
                                                          512.0 / 34816.0,
                                                          1024.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_2(std::make_tuple(189.0 / 34816.0,
                                                          210.0 / 34816.0,
                                                          240.0 / 34816.0,
                                                          288.0 / 34816.0,
                                                          284.0 / 34816.0,
                                                          768.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_3(std::make_tuple(175.0 / 34816.0,
                                                          200.0 / 34816.0,
                                                          240.0 / 34816.0,
                                                          320.0 / 34816.0,
                                                          640.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_4(std::make_tuple(175.0 / 34816.0,
                                                          210.0 / 34816.0,
                                                          280.0 / 34816.0,
                                                          560.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_5(std::make_tuple(189.0 / 34816.0,
                                                          252.0 / 34816.0,
                                                          504.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_6(std::make_tuple(231.0 / 34816.0,
                                                          462.0 / 34816.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_7_7(std::make_tuple(429.0 / 34816.0));

PolynomialInMonomialBasis<double, double, 8, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_0(std::make_tuple(6435.0 / 622592.0,
                                                          6864.0 / 622592.0,
                                                          7392.0 / 622592.0,
                                                          8064.0 / 622592.0,
                                                          8960.0 / 622592.0,
                                                          10240.0 / 622592.0,
                                                          12288.0 / 622592.0,
                                                          16384.0 / 622592.0,
                                                          32768.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_1(std::make_tuple(3432.0 / 622592.0,
                                                          3696.0 / 622592.0,
                                                          4032.0 / 622592.0,
                                                          4480.0 / 622592.0,
                                                          5120.0 / 622592.0,
                                                          6144.0 / 622592.0,
                                                          8192.0 / 622592.0,
                                                          16384.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_2(std::make_tuple(2772.0 / 622592.0,
                                                          3024.0 / 622592.0,
                                                          3360.0 / 622592.0,
                                                          3840.0 / 622592.0,
                                                          4608.0 / 622592.0,
                                                          6144.0 / 622592.0,
                                                          12288.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_3(std::make_tuple(2520.0 / 622592.0,
                                                          2800.0 / 622592.0,
                                                          3200.0 / 622592.0,
                                                          3840.0 / 622592.0,
                                                          5120.0 / 622592.0,
                                                          10240.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_4(std::make_tuple(2450.0 / 622592.0,
                                                          2800.0 / 622592.0,
                                                          3360.0 / 622592.0,
                                                          4480.0 / 622592.0,
                                                          8960.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_5(std::make_tuple(2520.0 / 622592.0,
                                                          3024.0 / 622592.0,
                                                          4032.0 / 622592.0,
                                                          8064.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_6(std::make_tuple(2772.0 / 622592.0,
                                                          3696.0 / 622592.0,
                                                          7392.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_7(std::make_tuple(3432.0 / 622592.0,
                                                          6864.0 / 622592.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_8_8(std::make_tuple(6435.0 / 622592.0));

PolynomialInMonomialBasis<double, double, 9, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_0(std::make_tuple(12155.0 / 1376256.0,
                                                          12870.0 / 1376256.0,
                                                          13728.0 / 1376256.0,
                                                          14784.0 / 1376256.0,
                                                          16128.0 / 1376256.0,
                                                          17920.0 / 1376256.0,
                                                          20480.0 / 1376256.0,
                                                          24576.0 / 1376256.0,
                                                          32768.0 / 1376256.0,
                                                          65536.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 8, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_1(std::make_tuple(6435.0 / 1376256.0,
                                                          6864.0 / 1376256.0,
                                                          7392.0 / 1376256.0,
                                                          8064.0 / 1376256.0,
                                                          8960.0 / 1376256.0,
                                                          10240.0 / 1376256.0,
                                                          12288.0 / 1376256.0,
                                                          16384.0 / 1376256.0,
                                                          32768.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_2(std::make_tuple(5148.0 / 1376256.0,
                                                          5544.0 / 1376256.0,
                                                          6048.0 / 1376256.0,
                                                          6720.0 / 1376256.0,
                                                          7680.0 / 1376256.0,
                                                          9216.0 / 1376256.0,
                                                          12288.0 / 1376256.0,
                                                          24576.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_3(std::make_tuple(4620.0 / 1376256.0,
                                                          5040.0 / 1376256.0,
                                                          5600.0 / 1376256.0,
                                                          6400.0 / 1376256.0,
                                                          7680.0 / 1376256.0,
                                                          10240.0 / 1376256.0,
                                                          20480.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_4(std::make_tuple(4410.0 / 1376256.0,
                                                          4900.0 / 1376256.0,
                                                          5600.0 / 1376256.0,
                                                          6720.0 / 1376256.0,
                                                          8960.0 / 1376256.0,
                                                          17920.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_5(std::make_tuple(4410.0 / 1376256.0,
                                                          5040.0 / 1376256.0,
                                                          6048.0 / 1376256.0,
                                                          8064.0 / 1376256.0,
                                                          16128.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_6(std::make_tuple(4620.0 / 1376256.0,
                                                          5544.0 / 1376256.0,
                                                          7392.0 / 1376256.0,
                                                          14784.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_7(std::make_tuple(5148.0 / 1376256.0,
                                                          6864.0 / 1376256.0,
                                                          13728.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_8(std::make_tuple(6435.0 / 1376256.0,
                                                          12870.0 / 1376256.0));
PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
    fukushima_elliptic_Js_maclaurin_n_9_9(std::make_tuple(12155.0 / 1376256.0));

// A generator for the Maclaurin series for Fukushima's T function.
template<int n, template<typename, typename, int> class Evaluator>
class FukushimaTMaclaurin {
  template<typename>
  struct Generator;

  template<int... k>
  struct Generator<std::index_sequence<k...>> {
    static auto constexpr series = std::make_tuple(1.0 / (2.0 * k + 1.0)...);
  };

 public:
  static inline PolynomialInMonomialBasis<double, double, n, Evaluator> const
      polynomial{Generator<std::make_index_sequence<n + 1>>::series};
};

using FukushimaTMaclaurin1 = FukushimaTMaclaurin<1, HornerEvaluator>;
using FukushimaTMaclaurin2 = FukushimaTMaclaurin<2, HornerEvaluator>;
using FukushimaTMaclaurin3 = FukushimaTMaclaurin<3, HornerEvaluator>;
using FukushimaTMaclaurin4 = FukushimaTMaclaurin<4, EstrinEvaluator>;
using FukushimaTMaclaurin5 = FukushimaTMaclaurin<5, EstrinEvaluator>;
using FukushimaTMaclaurin6 = FukushimaTMaclaurin<6, EstrinEvaluator>;
using FukushimaTMaclaurin7 = FukushimaTMaclaurin<7, EstrinEvaluator>;
using FukushimaTMaclaurin8 = FukushimaTMaclaurin<8, EstrinEvaluator>;
using FukushimaTMaclaurin9 = FukushimaTMaclaurin<9, EstrinEvaluator>;
using FukushimaTMaclaurin10 = FukushimaTMaclaurin<10, EstrinEvaluator>;
using FukushimaTMaclaurin11 = FukushimaTMaclaurin<11, EstrinEvaluator>;
using FukushimaTMaclaurin12 = FukushimaTMaclaurin<12, EstrinEvaluator>;

// Polynomials for EllipticK.  The last part of the name indicates the value of
// m around which the approximation is valid.
// TODO(phl): Some of these polynomials use the Horner scheme because the Estrin
// scheme introduces inaccuracies.  Investigate why.
PolynomialInMonomialBasis<double, double, 10, HornerEvaluator> const
    elliptic_K_taylor_0_05(std::make_tuple(1.591003453790792180,
                                           0.416000743991786912,
                                           0.245791514264103415,
                                           0.179481482914906162,
                                           0.144556057087555150,
                                           0.123200993312427711,
                                           0.108938811574293531,
                                           0.098853409871592910,
                                           0.091439629201749751,
                                           0.085842591595413900,
                                           0.081541118718303215));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    elliptic_K_taylor_0_15(std::make_tuple(1.635256732264579992,
                                           0.471190626148732291,
                                           0.309728410831499587,
                                           0.252208311773135699,
                                           0.226725623219684650,
                                           0.215774446729585976,
                                           0.213108771877348910,
                                           0.216029124605188282,
                                           0.223255831633057896,
                                           0.234180501294209925,
                                           0.248557682972264071,
                                           0.266363809892617521,
                                           0.287728452156114668));
PolynomialInMonomialBasis<double, double, 11, EstrinEvaluator> const
    elliptic_K_taylor_0_25(std::make_tuple(1.685750354812596043,
                                           0.541731848613280329,
                                           0.401524438390690257,
                                           0.369642473420889090,
                                           0.376060715354583645,
                                           0.405235887085125919,
                                           0.453294381753999079,
                                           0.520518947651184205,
                                           0.609426039204995055,
                                           0.724263522282908870,
                                           0.871013847709812357,
                                           1.057652872753547036));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    elliptic_K_taylor_0_35(std::make_tuple(1.744350597225613243,
                                           0.634864275371935304,
                                           0.539842564164445538,
                                           0.571892705193787391,
                                           0.670295136265406100,
                                           0.832586590010977199,
                                           1.073857448247933265,
                                           1.422091460675497751,
                                           1.920387183402304829,
                                           2.632552548331654201,
                                           3.652109747319039160,
                                           5.115867135558865806,
                                           7.224080007363877411));
PolynomialInMonomialBasis<double, double, 13, EstrinEvaluator> const
    elliptic_K_taylor_0_45(std::make_tuple(1.813883936816982644,
                                           0.763163245700557246,
                                           0.761928605321595831,
                                           0.951074653668427927,
                                           1.315180671703161215,
                                           1.928560693477410941,
                                           2.937509342531378755,
                                           4.594894405442878062,
                                           7.330071221881720772,
                                           11.87151259742530180,
                                           19.45851374822937738,
                                           32.20638657246426863,
                                           53.73749198700554656,
                                           90.27388602940998849));
PolynomialInMonomialBasis<double, double, 14, HornerEvaluator> const
    elliptic_K_taylor_0_55(std::make_tuple(1.898924910271553526,
                                           0.950521794618244435,
                                           1.151077589959015808,
                                           1.750239106986300540,
                                           2.952676812636875180,
                                           5.285800396121450889,
                                           9.832485716659979747,
                                           18.78714868327559562,
                                           36.61468615273698145,
                                           72.45292395127771801,
                                           145.1079577347069102,
                                           293.4786396308497026,
                                           598.3851815055010179,
                                           1228.420013075863451,
                                           2536.529755382764488));
PolynomialInMonomialBasis<double, double, 16, EstrinEvaluator> const
    elliptic_K_taylor_0_65(std::make_tuple(2.007598398424376302,
                                           1.248457231212347337,
                                           1.926234657076479729,
                                           3.751289640087587680,
                                           8.119944554932045802,
                                           18.66572130873555361,
                                           44.60392484291437063,
                                           109.5092054309498377,
                                           274.2779548232413480,
                                           697.5598008606326163,
                                           1795.716014500247129,
                                           4668.381716790389910,
                                           12235.76246813664335,
                                           32290.17809718320818,
                                           85713.07608195964685,
                                           228672.1890493117096,
                                           612757.2711915852774));
PolynomialInMonomialBasis<double, double, 19, EstrinEvaluator> const
    elliptic_K_taylor_0_75(std::make_tuple(2.156515647499643235,
                                           1.791805641849463243,
                                           3.826751287465713147,
                                           10.38672468363797208,
                                           31.40331405468070290,
                                           100.9237039498695416,
                                           337.3268282632272897,
                                           1158.707930567827917,
                                           4060.990742193632092,
                                           14454.00184034344795,
                                           52076.66107599404803,
                                           189493.6591462156887,
                                           695184.5762413896145,
                                           2.567994048255284686e6,
                                           9.541921966748386322e6,
                                           3.563492744218076174e7,
                                           1.336692984612040871e8,
                                           5.033521866866284541e8,
                                           1.901975729538660119e9,
                                           7.208915015330103756e9));
PolynomialInMonomialBasis<double, double, 15, EstrinEvaluator> const
    elliptic_K_taylor_0_825(std::make_tuple(2.318122621712510589,
                                            2.616920150291232841,
                                            7.897935075731355823,
                                            30.50239715446672327,
                                            131.4869365523528456,
                                            602.9847637356491617,
                                            2877.024617809972641,
                                            14110.51991915180325,
                                            70621.44088156540229,
                                            358977.2665825309926,
                                            1.847238263723971684e6,
                                            9.600515416049214109e6,
                                            5.030767708502366879e7,
                                            2.654441886527127967e8,
                                            1.408862325028702687e9,
                                            7.515687935373774627e9));
PolynomialInMonomialBasis<double, double, 19, EstrinEvaluator> const
    elliptic_K_taylor_0_875(std::make_tuple(2.473596173751343912,
                                            3.727624244118099310,
                                            15.60739303554930496,
                                            84.12850842805887747,
                                            506.9818197040613935,
                                            3252.277058145123644,
                                            21713.24241957434256,
                                            149037.0451890932766,
                                            1.043999331089990839e6,
                                            7.427974817042038995e6,
                                            5.350383967558661151e7,
                                            3.892498869948708474e8,
                                            2.855288351100810619e9,
                                            2.109007703876684053e10,
                                            1.566998339477902014e11,
                                            1.170222242422439893e12,
                                            8.777948323668937971e12,
                                            6.610124275248495041e13,
                                            4.994880537133887989e14,
                                            3.785974339724029920e15));

// Polynomials for FukushimaEllipticBD.  The last part of the name indicates the
// value of m around which the approximation is valid.
PolynomialInMonomialBasis<double, double, 11, EstrinEvaluator> const
    fukushima_B_taylor_0_05(
        std::make_tuple(0.790401413584395132310045630540381158921005,
                        0.102006266220019154892513446364386528537788,
                        0.039878395558551460860377468871167215878458,
                        0.021737136375982167333478696987134316809322,
                        0.013960979767622057852185340153691548520857,
                        0.009892518822669142478846083436285145400444,
                        0.007484612400663335676130416571517444936951,
                        0.005934625664295473695080715589652011420808,
                        0.004874249053581664096949448689997843978535,
                        0.004114606930310886136960940893002069423559,
                        0.003550452989196176932747744728766021440856,
                        0.003119229959988474753291950759202798352266));
PolynomialInMonomialBasis<double, double, 11, EstrinEvaluator> const
    fukushima_D_taylor_0_05(
        std::make_tuple(0.800602040206397047799296975176819811774784,
                        0.313994477771767756849615832867393028789057,
                        0.205913118705551954501930953451976374435088,
                        0.157744346538923994475225014971416837073598,
                        0.130595077319933091909091103101366509387938,
                        0.113308474489758568672985167742047066367053,
                        0.101454199173630195376251916342483192174927,
                        0.0929187842072974367037702927967784464949434,
                        0.0865653801481680871714054745336652101162894,
                        0.0817279846651030135350056216958053404884715,
                        0.0779906657291070378163237851392095284454654,
                        0.075080426851268007156477347905308063808848));
PolynomialInMonomialBasis<double, double, 11, EstrinEvaluator> const
    fukushima_B_taylor_0_15(
        std::make_tuple(0.80102406445284489393880821604009991524037,
                        0.11069534452963401497502459778015097487115,
                        0.047348746716993717753569559936346358937777,
                        0.028484367255041422845322166419447281776162,
                        0.020277811444003597057721308432225505126013,
                        0.015965005853099119442287313909177068173564,
                        0.013441320273553634762716845175446390822633,
                        0.011871565736951439501853534319081030547931,
                        0.010868363672485520630005005782151743785248,
                        0.010231587232710564565903812652581252337697,
                        0.009849585546666211201566987057592610884309,
                        0.009656606347153765129943681090056980586986));
PolynomialInMonomialBasis<double, double, 11, EstrinEvaluator> const
    fukushima_D_taylor_0_15(
        std::make_tuple(0.834232667811735098431315595374145207701720,
                        0.360495281619098275577215529302260739976126,
                        0.262379664114505869328637749459234348287432,
                        0.223723944518094276386520735054801578584350,
                        0.206447811775681052682922746753795148394463,
                        0.199809440876486856438050774316751253389944,
                        0.199667451603795274869211409350873244844882,
                        0.204157558868236842039815028663379643303565,
                        0.212387467960572375038025392458549025660994,
                        0.223948914061499360356873401571821627069173,
                        0.238708097425597860161720875806632864507536,
                        0.256707203545463755643710021815937785120030));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_B_taylor_0_25(
        std::make_tuple(0.81259777291992049322557009456643357559904,
                        0.12110961794551011284012693733241967660542,
                        0.057293376831239877456538980381277010644332,
                        0.038509451602167328057004166642521093142114,
                        0.030783430301775232744816612250699163538318,
                        0.027290564934732526869467118496664914274956,
                        0.025916369289445198731886546557337255438083,
                        0.025847203343361799141092472018796130324244,
                        0.026740923539348854616932735567182946385269,
                        0.028464314554825704963640157657034405579849,
                        0.030995446237278954096189768338119395563447,
                        0.034384369179940975864103666824736551261799,
                        0.038738002072493935952384233588242422046537));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_D_taylor_0_25(
        std::make_tuple(0.873152581892675549645633563232643413901757,
                        0.420622230667770215976919792378536040460605,
                        0.344231061559450379368201151870166692934830,
                        0.331133021818721761888662390999045979071436,
                        0.345277285052808411877017306497954757532251,
                        0.377945322150393391759797943135325823338761,
                        0.427378012464553880508348757311170776829930,
                        0.494671744307822405584118022550673740404732,
                        0.582685115665646200824237214098764913658889,
                        0.695799207728083164790111837174250683834359,
                        0.840018401472533403272555302636558338772258,
                        1.023268503573606060588689738498395211300483,
                        1.255859085136282496149035687741403985044122));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_B_taylor_0_35(
        std::make_tuple(0.8253235579835158949845697805395190063745,
                        0.1338621160836877898575391383950840569989,
                        0.0710112935979886745743770664203746758134,
                        0.0541784774173873762208472152701393154906,
                        0.0494517449481029932714386586401273353617,
                        0.0502221962241074764652127892365024315554,
                        0.0547429131718303528104722303305931350375,
                        0.0627462579270016992000788492778894700075,
                        0.0746698810434768864678760362745179321956,
                        0.0914808451777334717996463421986810092918,
                        0.1147050921109978235104185800057554574708,
                        0.1465711325814398757043492181099197917984,
                        0.1902571373338462844225085057953823854177));
PolynomialInMonomialBasis<double, double, 13, EstrinEvaluator> const
    fukushima_D_taylor_0_35(
        std::make_tuple(0.9190270392420973478848471774160778462738,
                        0.5010021592882475139767453081737767171354,
                        0.4688312705664568629356644841691659415972,
                        0.5177142277764000147059587510833317474467,
                        0.6208433913173031070711926900889045286988,
                        0.7823643937868697229213240489900179142670,
                        1.0191145350761029126165253557593691585239,
                        1.3593452027484960522212885423056424704073,
                        1.8457173023588279422916645725184952058635,
                        2.5410717031539207287662105618152273788399,
                        3.5374046552080413366422791595672470037341,
                        4.9692960029774259303491034652093672488707,
                        7.0338228700300311264031522795337352226926,
                        10.020043225034471401553194050933390974016));
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_B_taylor_0_45(
        std::make_tuple(0.8394795702706129706783934654948360410325,
                        0.1499164403063963359478614453083470750543,
                        0.0908319358194288345999005586556105610069,
                        0.0803470334833417864262134081954987019902,
                        0.0856384405004704542717663971835424473169,
                        0.1019547259329903716766105911448528069506,
                        0.1305748115336160150072309911623351523284,
                        0.1761050763588499277679704537732929242811,
                        0.2468351644029554468698889593583314853486,
                        0.3564244768677188553323196975301769697977,
                        0.5270025622301027434418321205779314762241,
                        0.7943896342593047502260866957039427731776,
                        1.2167625324297180206378753787253096783993));
PolynomialInMonomialBasis<double, double, 15, EstrinEvaluator> const
    fukushima_D_taylor_0_45(
        std::make_tuple(0.9744043665463696730314687662723484085813,
                        0.6132468053941609101234053415051402349752,
                        0.6710966695021669963502789954058993004082,
                        0.8707276201850861403618528872292437242726,
                        1.2295422312026907609906452348037196571302,
                        1.8266059675444205694817638548699906990301,
                        2.8069345309977627400322167438821024032409,
                        4.4187893290840281339827573139793805587268,
                        7.0832360574787653249799018590860687062869,
                        11.515088120557582942290563338274745712174,
                        18.931511185999274639516732819605594455165,
                        31.411996938204963878089048091424028309798,
                        52.520729454575828537934780076768577185134,
                        88.384854735065298062125622417251073520996,
                        149.56637449398047835236703116483092644714,
                        254.31790843104117434615624121937495622372));
PolynomialInMonomialBasis<double, double, 13, EstrinEvaluator> const
    fukushima_B_taylor_0_55(
        std::make_tuple(0.8554696151564199914087224774321783838373,
                        0.1708960726897395844132234165994754905373,
                        0.1213352290269482260207667564010437464156,
                        0.1282018835749474096272901529341076494573,
                        0.1646872814515275597348427294090563472179,
                        0.2374189087493817423375114793658754489958,
                        0.3692081047164954516884561039890508294508,
                        0.6056587338479277173311618664015401963868,
                        1.0337055615578127436826717513776452476106,
                        1.8189884893632678849599091011718520567105,
                        3.2793776512738509375806561547016925831128,
                        6.0298883807175363312261449542978750456611,
                        11.269796855577941715109155203721740735793,
                        21.354577850382834496786315532111529462693));
PolynomialInMonomialBasis<double, double, 16, EstrinEvaluator> const
    fukushima_D_taylor_0_55(
        std::make_tuple(1.04345529511513353426326823569160142342838,
                        0.77962572192850485048535711388072271372632,
                        1.02974236093206758187389128668777397528702,
                        1.62203722341135313022433907993860147395972,
                        2.78798953118534762046989770119382209443756,
                        5.04838148737206914685643655935236541332892,
                        9.46327761194348429539987572314952029503864,
                        18.1814899494276679043749394081463811247757,
                        35.5809805911791687037085198750213045708148,
                        70.6339354619144501276254906239838074917358,
                        141.828580083433059297030133195739832297859,
                        287.448751250132166257642182637978103762677,
                        587.115384649923076181773192202238389711345,
                        1207.06543522548061603657141890778290399603,
                        2495.58872724866422273012188618178997342537,
                        5184.69242939480644062471334944523925163600,
                        10817.2133369041327524988910635205356016939));
PolynomialInMonomialBasis<double, double, 15, EstrinEvaluator> const
    fukushima_B_taylor_0_65(
        std::make_tuple(0.8739200618486431359820482173294324246058,
                        0.1998140574823769459497418213885348159654,
                        0.1727696158780152128147094051876565603862,
                        0.2281069132842021671319791750725846795701,
                        0.3704681411180712197627619157146806221767,
                        0.6792712528848205545443855883980014994450,
                        1.3480084966817573020596179874311042267679,
                        2.8276709768538207038046918250872679902352,
                        6.1794682501239140840906583219887062092430,
                        13.935686010342811497608625663457407447757,
                        32.218929281059722026322932181420383764028,
                        76.006962959226101026399085304912635262362,
                        182.32144908775406957609058046006949657416,
                        443.51507644112648158679360783118806161062,
                        1091.8547229028388292980623647414961662223,
                        2715.7658664038195881056269799613407111521));
PolynomialInMonomialBasis<double, double, 17, EstrinEvaluator> const
    fukushima_D_taylor_0_65(
        std::make_tuple(1.13367833657573316571671258513452768536080,
                        1.04864317372997039116746991765351150490010,
                        1.75346504119846451588826580872136305225406,
                        3.52318272680338551269021618722443199230946,
                        7.74947641381397458240336356601913534598302,
                        17.9864500558507330560532617743406294626849,
                        43.2559163462326133313977294448984936591235,
                        106.681534454096017031613223924991564429656,
                        268.098486573117433951562111736259672695883,
                        683.624114850289804796762005964155730439745,
                        1763.49708521918740723028849567007874329637,
                        4592.37475383116380899419201719007475759114,
                        12053.4410190488892782190764838488156555734,
                        31846.6630207420816960681624497373078887317,
                        84621.2213590568080177035346867495326879117,
                        225956.423182907889987641304430180593010940,
                        605941.517281758859958050194535269219533685,
                        1.63108259953926832083633544697688841456604e6));
PolynomialInMonomialBasis<double, double, 18, EstrinEvaluator> const
    fukushima_B_taylor_0_75(
        std::make_tuple(0.895902820924731621258525533131864225704,
                        0.243140003766786661947749288357729051637,
                        0.273081875594105531575351304277604081620,
                        0.486280007533573323895498576715458103274,
                        1.082747437228230914750752674136983406683,
                        2.743445290986452500459431536349945437824,
                        7.555817828670234627048618342026400847824,
                        22.05194082493752427472777448620986154515,
                        67.15640644740229407624192175802742979626,
                        211.2722537881770961691291434845898538537,
                        681.9037843053270682273212958093073895805,
                        2246.956231592536516768812462150619631201,
                        7531.483865999711792004783423815426725079,
                        25608.51260130241579018675054866136922157,
                        88140.74740089604971425934283371277143256,
                        306564.4242098446591430938434419151070722,
                        1.076036077811072193752770590363885180738e6,
                        3.807218502573632648224286313875985190526e6,
                        1.356638224422139551020110323739879481042e7));
PolynomialInMonomialBasis<double, double, 20, EstrinEvaluator> const
    fukushima_D_taylor_0_75(
        std::make_tuple(1.26061282657491161418014946566845780315983,
                        1.54866563808267658056930177790599939977154,
                        3.55366941187160761540650011660758187283401,
                        9.90044467610439875577300608183010716301714,
                        30.3205666174524719862025105898574414438275,
                        98.1802586588830891484913119780870074464833,
                        329.771010434557055036273670551546757245808,
                        1136.65598974289039303581967838947708073239,
                        3993.83433574622979757935610692842933356144,
                        14242.7295865552708506232731633468180669284,
                        51394.7572916887209594591528374806790960057,
                        187246.702914623152141768788258141932569037,
                        687653.092375389902708761221294282367947659,
                        2.54238553565398227033448846432182516906624e6,
                        9.45378121934749027243313241962076028066811e6,
                        3.53283630179709170835024033154326126569613e7,
                        1.32593262383393014923560730485845833322771e8,
                        4.99544968184054821463279808395426941549833e8,
                        1.88840934729443872364972817525484292678543e9,
                        7.16026753447893719179055010636502508063102e9,
                        2.72233079469633962247554894093591262281929e10));
PolynomialInMonomialBasis<double, double, 14, EstrinEvaluator> const
    fukushima_B_taylor_0_825(
        std::make_tuple(0.915922052601931494319853880201442948834592,
                        0.294714252429483394379515488141632749820347,
                        0.435776709264636140422971598963772380161131,
                        1.067328246493644238508159085364429570207744,
                        3.327844118563268085074646976514979307993733,
                        11.90406004445092906188837729711173326621810,
                        46.47838820224626393512400481776284680677096,
                        192.7556002578809476962739389101964074608802,
                        835.3356299261900063712302517586717381557137,
                        3743.124548343029102644419963712353854902019,
                        17219.07731004063094108708549153310467326395,
                        80904.60401669850158353080543152212152282878,
                        386808.3292751742460123683674607895217760313,
                        1.876487670110449342170327796786290400635732e6,
                        9.216559908641567755240142886998737950775910e6));
PolynomialInMonomialBasis<double, double, 17, EstrinEvaluator> const
    fukushima_D_taylor_0_825(
        std::make_tuple(1.402200569110579095046054435635136986038164,
                        2.322205897861749446534352741005347103992773,
                        7.462158366466719682730245467372788273333992,
                        29.43506890797307903104978364254987042421285,
                        128.1590924337895775262509354898066132182429,
                        591.0807036911982326384997979640812493154316,
                        2830.546229607726377048576057043685514661188,
                        13917.76431889392229954434840686741305556862,
                        69786.10525163921228258055074102587429394212,
                        355234.1420341879634781808899208309503519936,
                        1.830019186413931053503912913904321703777885e6,
                        9.519610812032515607466102200648641326190483e6,
                        4.992086875574849453986274042758566713803723e7,
                        2.635677009826023473846461512029006874800883e8,
                        1.399645765120061118824228996253541612110338e9,
                        7.469935792837635004663183580452618726280406e9,
                        4.004155595835610574316003488168804738481448e10,
                        2.154630668144966654449602981243932210695662e11));
PolynomialInMonomialBasis<double, double, 18, EstrinEvaluator> const
    fukushima_B_taylor_0_875(
        std::make_tuple(0.931906061029524827613331428871579482766771,
                        0.348448029538453860999386797137074571589376,
                        0.666809178846938247558793864839434184202736,
                        2.210769135708128662563678717558631573758222,
                        9.491765048913406881414290930355300611703187,
                        47.09304791027740853381457907791343619298913,
                        255.9200460211233087050940506395442544885608,
                        1480.029532675805407554800779436693505109703,
                        8954.040904734313578374783155553041065984547,
                        56052.48220982686949967604699243627330816542,
                        360395.7241626000916973524840479780937869149,
                        2.367539415273216077520928806581689330885103e6,
                        1.582994957277684102454906900025484391190264e7,
                        1.074158093278511100137056972128875270067228e8,
                        7.380585460239595691878086073095523043390649e8,
                        5.126022002555101496684687154904781856830296e9,
                        3.593534065502416588712409180013118409428367e10,
                        2.539881257612812212004146637239987308133582e11,
                        1.808180007145359569674767150594344316702507e12));
PolynomialInMonomialBasis<double, double, 20, EstrinEvaluator> const
    fukushima_D_taylor_0_875(
        std::make_tuple(1.541690112721819084362258323861459983048179,
                        3.379176214579645449453938918349243359477706,
                        14.94058385670236671625328259137998668324435,
                        81.91773929235074880784578753539752529822986,
                        497.4900546551479866036061853049402721939835,
                        3205.184010234846235275447901572262470252768,
                        21457.32237355321925571253220641357074594515,
                        147557.0156564174712105689758692929775004292,
                        1.035045290185256525452269053775273002725343e6,
                        7.371922334832212125197513363695905834126154e6,
                        5.314344395142401141792228169170505958906345e7,
                        3.868823475795976312985118115567305767603128e8,
                        2.839458401528033778425531336599562337200510e9,
                        2.098266122943898941547136470383199468548861e10,
                        1.559617754017662417944194874282275405676282e11,
                        1.165096220419884791236699872205721392201682e12,
                        8.742012983013913804987431275193291316808766e12,
                        6.584725462672366918676967847406180155459650e13,
                        4.976798737062434393396993620379481464465749e14,
                        3.773018634056605404718444239040628892506293e15,
                        2.868263194837819660109735981973458220407767e16));

// NOTE(phl): The following polynomials differ slightly from the original code
// but they match more closely those in [Fuk11a].  The notation follows
// [Fuk11a].
// A polynomial for B٭X(m) / m.
PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator> const
    fukushima_B٭X_maclaurin(std::make_tuple(-1.0 / 4.0,
                                            -1.0 / 32.0,
                                            -3.0 / 256.0,
                                            -25.0 / 4096.0,
                                            -245.0 / 65536.0,
                                            -1323.0 / 524288.0,
                                            -7623.0 / 4194304.0,
                                            -184041.0 / 134217728.0));

// A polynomial for EX(m) / m.
PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator> const
    fukushima_EX_maclaurin(std::make_tuple(1.0 / 4.0,
                                           3.0 / 32.0,
                                           15.0 / 256.0,
                                           175.0 / 4096.0,
                                           2205.0 / 65536.0,
                                           14553.0 / 524288.0,
                                           99099.0 / 4194304.0,
                                           2760615.0 / 134217728.0));

// A polynomial for KX(m).
PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator> const
    fukushima_KX_maclaurin(std::make_tuple(1.0 / 2.0,
                                           1.0 / 8.0,
                                           9.0 / 128.0,
                                           25.0 / 512.0,
                                           1225.0 / 32768.0,
                                           3969.0 / 131072.0,
                                           53361.0 / 2097152.0,
                                           184041.0 / 8388608.0));

// A polynomial for EX(m).
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_EX_taylor_0_05(
        std::make_tuple(0.02548395442966088473597712420249483947953 * 0.5,
                        0.51967384324140471318255255900132590084179 * 0.5,
                        0.20644951110163173131719312525729037023377 * 0.5,
                        0.13610952125712137420240739057403788152260 * 0.5,
                        0.10458014040566978574883406877392984277718 * 0.5,
                        0.08674612915759188982465635633597382093113 * 0.5,
                        0.07536380269617058326770965489534014190391 * 0.5,
                        0.06754544594618781950496091910264174396541 * 0.5,
                        0.06190939688096410201497509102047998554900 * 0.5,
                        0.05771071515451786553160533778648705873199 * 0.5,
                        0.05451217098672207169493767625617704078257 * 0.5,
                        0.05204028407582600387265992107877094920787 * 0.5,
                        0.05011532514520838441892567405879742720039 * 0.5));

// A polynomial for KX(m).
PolynomialInMonomialBasis<double, double, 12, EstrinEvaluator> const
    fukushima_KX_taylor_0_05(std::make_tuple(
        (1.0 + 0.01286425658832983978282698630501405107893) * 0.5,
        0.26483429894479586582278131697637750604652 * 0.5,
        0.15647573786069663900214275050014481397750 * 0.5,
        0.11426146079748350067910196981167739749361 * 0.5,
        0.09202724415743445309239690377424239940545 * 0.5,
        0.07843218831801764082998285878311322932444 * 0.5,
        0.06935260142642158347117402021639363379689 * 0.5,
        0.06293203529021269706312943517695310879457 * 0.5,
        0.05821227592779397036582491084172892108196 * 0.5,
        0.05464909112091564816652510649708377642504 * 0.5,
        0.05191068843704411873477650167894906357568 * 0.5,
        0.04978344771840508342564702588639140680363 * 0.5,
        0.04812375496807025605361215168677991360500 * 0.5));

Angle BulirschCel(double kc, double const nc, double a, double b) {
  // These values should give us 14 digits of accuracy, see [Bul69].
  constexpr double ca = 1.0e-7;
  constexpr double kc_nearly_0 = 1.0e-14;

  // The identifiers below follow exactly [Bul69].  Note the (uncommon) use of
  // non-const parameters to mimic [Bul69].
  double p = nc;
  if (kc == 0.0) {
    if (b == 0.0) {
      kc = kc_nearly_0;
    } else {
      // "If in this case b ≠ 0 then cel is undefined."
      DLOG(ERROR) << "kc = " << kc << " nc = " << nc << " a = " << a
                  << " b = " << b;
      return std::numeric_limits<Angle>::quiet_NaN();
    }
  }
  kc = Abs(kc);
  double e = kc;
  double m = 1.0;

  // Initial values for p, a, b.
  if (p > 0.0) {
    p = Sqrt(p);
    b = b / p;
  } else {
    double f = kc * kc;
    double q = 1.0 - f;
    double g = 1.0 - p;
    f = f - p;
    q = (b - a * p) * q;
    p = Sqrt(f / g);
    a = (a - b) / g;
    b = a * p - q / (g * g * p);
  }

  // Bartky's algorithm.
  for (;;) {
    double f = a;
    a = b / p + a;
    double g = e / p;
    b = f * g + b;
    b = b + b;
    p = g + p;
    g = m;
    m = kc + m;
    if (Abs(g - kc) <= g * ca) {
      break;
    }
    kc = Sqrt(e);
    kc = kc + kc;
    e = kc * m;
  }
  return (π / 2) * (a * m + b) / (m * (m + p)) * Radian;
}

template<int degree>
double EllipticNomeQ(double const mc) {
  return mc *
         EllipticNomeQMaclaurin<degree - 1, EstrinEvaluator>::polynomial.
             Evaluate(mc);
}

void FukushimaEllipticBD(double const mc, Angle& B_m, Angle& D_m) {
  double const m = 1.0 - mc;
  if (m < std::numeric_limits<double>::epsilon() / 2.0) {
    B_m = π / 4 * Radian;
    D_m = π / 4 * Radian;
  } else if (mc < std::numeric_limits<double>::epsilon() / 2.0) {
    B_m = 1.0 * Radian;
    D_m = ((2.0 * log_2 - 1.0) - 0.5 * std::log(mc)) * Radian;
  } else if (mc < 0.1) {
    // This algorithm (from [Fuk18]) differs from the one in [Fuk11a] because it
    // divides log(q(mc)), not just log(mc / 16).  It tries to retain the same
    // notation, though.  See the documentation (Fukushima.pdf).
    double const nome = EllipticNomeQ<16>(mc);
    double const X_mc = -std::log(nome);  // X(mc).
    double KX_mc;  // KX(mc).
    double EX_mc;  // EX(mc).
    if (mc < 0.01) {
      KX_mc = fukushima_KX_maclaurin.Evaluate(mc);
      EX_mc = mc * fukushima_EX_maclaurin.Evaluate(mc);
    } else {
      double const mx = mc - 0.05;
      KX_mc = fukushima_KX_taylor_0_05.Evaluate(mx);
      EX_mc = fukushima_EX_taylor_0_05.Evaluate(mx);
    }
    // Equivalent to Fukushima’s code [Fuk18], but much simplified.
    double const one_over_two_KX_mc = 0.5 / KX_mc;
    B_m = (X_mc * (EX_mc - mc * KX_mc) + one_over_two_KX_mc) * Radian / m;
    D_m = X_mc * KX_mc * Radian - B_m;
  } else if (m <= 0.01) {
    B_m = (-π * Radian) * fukushima_B٭X_maclaurin.Evaluate(m);
    D_m = (π * Radian) * fukushima_EX_maclaurin.Evaluate(m);
  } else if (m <= 0.1) {
    double const mx = 0.95 - mc;
    B_m = fukushima_B_taylor_0_05.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_05.Evaluate(mx) * Radian;
  } else if (m <= 0.2) {
    double const mx = 0.85 - mc;
    B_m = fukushima_B_taylor_0_15.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_15.Evaluate(mx) * Radian;
  } else if (m <= 0.3) {
    double const mx = 0.75 - mc;
    B_m = fukushima_B_taylor_0_25.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_25.Evaluate(mx) * Radian;
  } else if (m <= 0.4) {
    double const mx = 0.65 - mc;
    B_m = fukushima_B_taylor_0_35.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_35.Evaluate(mx) * Radian;
  } else if (m <= 0.5) {
    double const mx = 0.55 - mc;
    B_m = fukushima_B_taylor_0_45.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_45.Evaluate(mx) * Radian;
  } else if (m <= 0.6) {
    double const mx = 0.45 - mc;
    B_m = fukushima_B_taylor_0_55.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_55.Evaluate(mx) * Radian;
  } else if (m <= 0.7) {
    double const mx = 0.35 - mc;
    B_m = fukushima_B_taylor_0_65.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_65.Evaluate(mx) * Radian;
  } else if (m <= 0.8) {
    double const mx = 0.25 - mc;
    B_m = fukushima_B_taylor_0_75.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_75.Evaluate(mx) * Radian;
  } else if (m <= 0.85) {
    double const mx = 0.175 - mc;
    B_m = fukushima_B_taylor_0_825.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_825.Evaluate(mx) * Radian;
  } else {
    double const mx = 0.125 - mc;
    B_m = fukushima_B_taylor_0_875.Evaluate(mx) * Radian;
    D_m = fukushima_D_taylor_0_875.Evaluate(mx) * Radian;
  }
}

template<typename ThirdKind, typename>
void FukushimaEllipticBDJ(double const nc,
                          double const mc,
                          Angle& B_m,
                          Angle& D_m,
                          ThirdKind& J_n_m) {
  if (mc < 0) {  // m > 1
    DLOG(FATAL) << "NYI";
  } else if (mc > 1) {  // m < 0
    double const m = 1 - mc;
    double const mN = -m / mc;
    double const mcN = 1 / mc;

    Angle B_mN{uninitialized};
    Angle D_mN{uninitialized};
    FukushimaEllipticBD(mcN, B_mN, D_mN);

    double const sqrt_mcN = Sqrt(mcN);
    B_m = sqrt_mcN * D_mN;
    D_m = sqrt_mcN * B_mN;
  } else {
    FukushimaEllipticBD(mc, B_m, D_m);
  }

  // According to [Fuk12], A.1, the Bulirsch function will work for all values
  // of m and n.
  if constexpr (should_compute<ThirdKind>) {
    // See [Bul69], special examples after equation (1.2.2).
    double const kc = Sqrt(mc);
    J_n_m = BulirschCel(kc, nc, /*a=*/0.0, /*b=*/1.0);
  }
}

// Note that the identifiers in the function definition are not the same as
// those in the function declaration.
// The declaration follows [Fuk11b], equations (9) and (10), and [Fuk12],
// equation (24), whereas the definition follows [Fuk11b], section 2.2, and
// [Fuk12], section 3.3.
// [Fuk11b] (for B and D) calls the index j and [Fuk12] (for J) calls it i; we
// use i everywhere.
template<typename ThirdKind, typename>
void FukushimaEllipticBcDcJc(double const c₀,
                             double const n,
                             double const mc,
                             Angle& Bᵢ,
                             Angle& Dᵢ,
                             ThirdKind& Jᵢ) {
  // See [Fuk11b] section 2.2 for the determination of xS.
  constexpr double xS = 0.1;
  // The maximum number of iterations in the first loop below.
  // NOTE(phl): I couldn't find a justification for this number.
  constexpr int max_transformations = 10;

  double y[max_transformations + 1];
  double s[max_transformations + 1];
  double cd[max_transformations + 1];

  double const m = 1.0 - mc;
  double const h = n * (1.0 - n) * (n - m);
  double const x₀ = c₀ * c₀;
  double const y₀ = 1.0 - x₀;

  // Alternate half and double argument transformations, when cancellations
  // would occur, [Fuk12] section 3.3.

  // Half argument transformation of c.
  y[0] = y₀;
  s[0] = Sqrt(y₀);
  double cᵢ = c₀;
  double xᵢ = x₀;
  int i = 0;  // Note that this variable is used after the loop.
  for (; xᵢ <= xS; ++i) {
    DCHECK_LT(i, max_transformations)
        << "c₀ = " << c₀ << " n = " << n << " mc = " << mc;
    double const dᵢ = Sqrt(mc + m * xᵢ);
    xᵢ = (cᵢ + dᵢ) / (1.0 + dᵢ);
    double const yᵢ = 1.0 - xᵢ;
    y[i + 1] = yᵢ;
    s[i + 1] = Sqrt(yᵢ);
    cd[i] = cᵢ * dᵢ;
    cᵢ = Sqrt(xᵢ);
  }
  int const I = i;  // The index at termination.

  // Switch to the normal algorithm.
  FukushimaEllipticBsDsJs(s[I], n, mc, Bᵢ, Dᵢ, Jᵢ);

  // Double argument transformation of B, D, J, [Fuk11b] equation (16) and
  // [Fuk12] equations (33–35).
  for (int i = I; i > 0; --i) {
    double const sy = s[i - 1] * y[i];
    Bᵢ = 2.0 * Bᵢ - sy * Radian;
    Dᵢ = Dᵢ + (Dᵢ + sy * Radian);
    if constexpr (should_compute<ThirdKind>) {
      double const t = sy / (1.0 - n * (y[i - 1] - y[i] * cd[i - 1]));
      Jᵢ = Jᵢ + (Jᵢ + FukushimaT(t, h));
    }
  }
}

// Note that the identifiers in the function definition are not the same as
// those in the function declaration.
// The declaration follows [Fuk11b], equations (9) and (10), and [Fuk12],
// equation (24), whereas the definition follows [Fuk11b], section 2.2, and
// [Fuk12], section 3.3.
// [Fuk11b] (for B and D) calls the index j and [Fuk12] (for J) calls it i; we
// use i everywhere.
template<typename ThirdKind, typename>
void FukushimaEllipticBsDsJs(double const s₀,
                             double const n,
                             double const mc,
                             Angle& Bᵢ,
                             Angle& Dᵢ,
                             ThirdKind& Jᵢ) {
  // See [Fuk12] section 3.5 for the determination of yB.
  constexpr double yB = 0.01622;
  // The maximum number of argument transformations, related to yB.  This is the
  // maximum number of iterations in the first loop below.
  constexpr int max_transformations = 10;

  double y[max_transformations + 1];
  double s[max_transformations + 1];
  double cd[max_transformations + 1];

  // Half and double argument transformations, [Fuk12] section 3.3.
  double const m = 1.0 - mc;
  double const h = n * (1.0 - n) * (n - m);
  double const y₀ = s₀ * s₀;

  // Half argument transformation of s.
  y[0] = y₀;
  s[0] = s₀;
  double yᵢ = y₀;
  int i = 0;  // Note that this variable is used after the loop.
  for (; yᵢ >= yB; ++i) {
    DCHECK_LT(i, max_transformations)
        << "s₀ = " << s₀ << " n = " << n << " mc = " << mc;
    double const cᵢ = Sqrt(1.0 - yᵢ);
    double const dᵢ = Sqrt(1.0 - m * yᵢ);
    yᵢ = yᵢ / ((1.0 + cᵢ) * (1.0 + dᵢ));
    y[i + 1] = yᵢ;
    s[i + 1] = Sqrt(yᵢ);
    cd[i] = cᵢ * dᵢ;
  }
  int const I = i;  // The index at termination.

  // Maclaurin series, [Fuk11b] equation (15) and [Fuk12] equation (32).
  Angle Σ_Bₗ_m_yˡ{uninitialized};
  Angle Σ_Dₗ_m_yˡ{uninitialized};
  FukushimaEllipticBsDsMaclaurinSeries(yᵢ, m, Σ_Bₗ_m_yˡ, Σ_Dₗ_m_yˡ);
  Bᵢ = s[I] * Σ_Bₗ_m_yˡ;
  Dᵢ = s[I] * yᵢ * Σ_Dₗ_m_yˡ;
  if constexpr (should_compute<ThirdKind>) {
    Jᵢ = s[I] * FukushimaEllipticJsMaclaurinSeries(yᵢ, n, m);
  }

  // Double argument transformation of B, D, J, [Fuk11b] equation (16) and
  // [Fuk12] equations (33–35).
  for (int i = I; i > 0; --i) {
    double const sy = s[i - 1] * y[i];
    Bᵢ = 2.0 * Bᵢ - sy * Radian;
    Dᵢ = Dᵢ + (Dᵢ + sy * Radian);
    if constexpr (should_compute<ThirdKind>) {
      double const t = sy / (1.0 - n * (y[i - 1] - y[i] * cd[i - 1]));
      Jᵢ = Jᵢ + (Jᵢ + FukushimaT(t, h));
    }
  }
}

// See [Fuk11b], section 2.3.
void FukushimaEllipticBsDsMaclaurinSeries(double const y,
                                          double const m,
                                          Angle& Σ_Bₗ_m_yˡ,
                                          Angle& Σ_Dₗ_m_yˡ) {
  double const Fs1 = FukushimaEllipticFsMaclaurin1::polynomial.Evaluate(m);
  double const Fs2 = FukushimaEllipticFsMaclaurin2::polynomial.Evaluate(m);
  double const Fs3 = FukushimaEllipticFsMaclaurin3::polynomial.Evaluate(m);
  double const Fs4 = FukushimaEllipticFsMaclaurin4::polynomial.Evaluate(m);
  double const Fs5 = FukushimaEllipticFsMaclaurin5::polynomial.Evaluate(m);
  double const Fs6 = FukushimaEllipticFsMaclaurin6::polynomial.Evaluate(m);
  double const Fs7 = FukushimaEllipticFsMaclaurin7::polynomial.Evaluate(m);
  double const Fs8 = FukushimaEllipticFsMaclaurin8::polynomial.Evaluate(m);
  double const Fs9 = FukushimaEllipticFsMaclaurin9::polynomial.Evaluate(m);
  double const Fs10 = FukushimaEllipticFsMaclaurin10::polynomial.Evaluate(m);
  double const Fs11 = FukushimaEllipticFsMaclaurin11::polynomial.Evaluate(m);

  auto const fukushima_elliptic_Ds_maclaurin =
      FukushimaEllipticDsBsMaclaurin<EstrinEvaluator>::MakeDsPolynomial(
      1.0, Fs1, Fs2, Fs3, Fs4, Fs5, Fs6, Fs7, Fs8, Fs9, Fs10, Fs11);
  Σ_Dₗ_m_yˡ = fukushima_elliptic_Ds_maclaurin.Evaluate(y) * Radian;

  auto const fukushima_elliptic_Bs_maclaurin =
      FukushimaEllipticDsBsMaclaurin<EstrinEvaluator>::MakeBsPolynomial(
      1.0, Fs1, Fs2, Fs3, Fs4, Fs5, Fs6, Fs7, Fs8, Fs9, Fs10, Fs11);
  Σ_Bₗ_m_yˡ = fukushima_elliptic_Bs_maclaurin.Evaluate(y) * Radian;
}

// See [Fuk12], section 3.4 and 3.5.
Angle FukushimaEllipticJsMaclaurinSeries(double const y,
                                         double const n,
                                         double const m) {
  // Maclaurin series in m whose coefficients are polynomials in n.  The index
  // is the degree in m (k in Fukushima's notation).
  PolynomialInMonomialBasis<double, double, 0, HornerEvaluator>
      fukushima_elliptic_Js_maclaurin_m_0(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_0_0.Evaluate(n)));
  PolynomialInMonomialBasis<double, double, 1, HornerEvaluator>
      fukushima_elliptic_Js_maclaurin_m_1(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_1_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_1_1.Evaluate(n)));
  PolynomialInMonomialBasis<double, double, 2, HornerEvaluator>
      fukushima_elliptic_Js_maclaurin_m_2(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_2_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_2_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_2_2.Evaluate(n)));
  PolynomialInMonomialBasis<double, double, 3, HornerEvaluator>
      fukushima_elliptic_Js_maclaurin_m_3(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_3_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_3_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_3_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_3_3.Evaluate(n)));
  PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_4(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_4_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_4_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_4_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_4_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_4_4.Evaluate(n)));
  // The first five coefficients Jₗ(n|m); the others are computed as needed.
  // TODO(egg): J₁(n|m) = 1/3; could we make it constexpr?
  double const J₁ = fukushima_elliptic_Js_maclaurin_m_0.Evaluate(m);
  double const J₂ = fukushima_elliptic_Js_maclaurin_m_1.Evaluate(m);
  double const J₃ = fukushima_elliptic_Js_maclaurin_m_2.Evaluate(m);
  double const J₄ = fukushima_elliptic_Js_maclaurin_m_3.Evaluate(m);
  double const J₅ = fukushima_elliptic_Js_maclaurin_m_4.Evaluate(m);
  // This function computes ∑ Jₗ(n|m) yˡ.  Since it has no constant term, we
  // compute it as y ∑ Jₗ(n|m) yˡ⁻¹.  In the remainder of this function,
  // |fukushima_elliptic_Js_maclaurin_y_*| is the polynomial ∑ Jₗ(n|m) yˡ⁻¹.
  if (y <= 6.0369310e-04) {
    // A Maclaurin series in y whose coefficients are polynomials in n and m.
    // The index is the degree in y of the series.  Since Js has no constant
    // term, this is (l - 1) in Fukushima's notation.
    PolynomialInMonomialBasis<double, double, 4, EstrinEvaluator>
        fukushima_elliptic_Js_maclaurin_y_4(
            std::make_tuple(J₁, J₂, J₃, J₄, J₅));
    return y * fukushima_elliptic_Js_maclaurin_y_4.Evaluate(y) * Radian;
  }

  PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_5(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_5_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_5_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_5_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_5_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_5_4.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_5_5.Evaluate(n)));
  double const J₆ = fukushima_elliptic_Js_maclaurin_m_5.Evaluate(m);
  if (y <= 2.0727505e-03) {
    PolynomialInMonomialBasis<double, double, 5, EstrinEvaluator>
        fukushima_elliptic_js_maclaurin_y_5(
            std::make_tuple(J₁, J₂, J₃, J₄, J₅, J₆));
    return y * fukushima_elliptic_js_maclaurin_y_5.Evaluate(y) * Radian;
  }

  PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_6(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_6_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_4.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_5.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_6_6.Evaluate(n)));
  double const J₇ = fukushima_elliptic_Js_maclaurin_m_6.Evaluate(m);
  if (y <= 5.0047026e-03) {
    PolynomialInMonomialBasis<double, double, 6, EstrinEvaluator>
        fukushima_elliptic_js_maclaurin_y_6(
            std::make_tuple(J₁, J₂, J₃, J₄, J₅, J₆, J₇));
    return y * fukushima_elliptic_js_maclaurin_y_6.Evaluate(y) * Radian;
  }

  PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_7(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_7_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_4.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_5.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_6.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_7_7.Evaluate(n)));
  double const J₈ = fukushima_elliptic_Js_maclaurin_m_7.Evaluate(m);
  if (y <= 9.6961652e-03) {
    PolynomialInMonomialBasis<double, double, 7, EstrinEvaluator>
        fukushima_elliptic_js_maclaurin_y_7(
            std::make_tuple(J₁, J₂, J₃, J₄, J₅, J₆, J₇, J₈));
    return y * fukushima_elliptic_js_maclaurin_y_7.Evaluate(y) * Radian;
  }

  PolynomialInMonomialBasis<double, double, 8, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_8(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_8_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_4.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_5.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_6.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_7.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_8_8.Evaluate(n)));
  double const J₉ = fukushima_elliptic_Js_maclaurin_m_8.Evaluate(m);
  if (y <= 1.6220210e-02) {
    PolynomialInMonomialBasis<double, double, 8, EstrinEvaluator>
        fukushima_elliptic_Js_maclaurin_y_8(
            std::make_tuple(J₁, J₂, J₃, J₄, J₅, J₆, J₇, J₈, J₉));
    return y * fukushima_elliptic_Js_maclaurin_y_8.Evaluate(y) * Radian;
  }

  PolynomialInMonomialBasis<double, double, 9, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_m_9(
          std::make_tuple(fukushima_elliptic_Js_maclaurin_n_9_0.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_1.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_2.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_3.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_4.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_5.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_6.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_7.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_8.Evaluate(n),
                          fukushima_elliptic_Js_maclaurin_n_9_9.Evaluate(n)));
  double const J₁₀ = fukushima_elliptic_Js_maclaurin_m_9.Evaluate(m);
  PolynomialInMonomialBasis<double, double, 9, EstrinEvaluator>
      fukushima_elliptic_Js_maclaurin_y_9(
          std::make_tuple(J₁, J₂, J₃, J₄, J₅, J₆, J₇, J₈, J₉, J₁₀));
  return y * fukushima_elliptic_Js_maclaurin_y_9.Evaluate(y) * Radian;
}

Angle FukushimaT(double const t, double const h) {
  double const z = -h * t * t;
  double const abs_z = abs(z);

  // NOTE(phl): One might be tempted to rewrite this statement using a binary
  // split of the interval [0, 1], but according to Table 1 of [Fuk12] the
  // distribution of z is very biased towards the small values, so this is
  // simpler and probably better.  (It also explains the position of z < 0 in
  // the list.)
  if (abs_z < 3.3306691e-16) {
    return t * Radian;
  } else if (abs_z < 2.3560805e-08) {
    return t * FukushimaTMaclaurin1::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 9.1939631e-06) {
    return t * FukushimaTMaclaurin2::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 1.7779240e-04) {
    return t * FukushimaTMaclaurin3::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 1.0407839e-03) {
    return t * FukushimaTMaclaurin4::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 3.3616998e-03) {
    return t * FukushimaTMaclaurin5::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 7.7408014e-03) {
    return t * FukushimaTMaclaurin6::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 1.4437181e-02) {
    return t * FukushimaTMaclaurin7::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 2.3407312e-02) {
    return t * FukushimaTMaclaurin8::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 3.4416203e-02) {
    return t * FukushimaTMaclaurin9::polynomial.Evaluate(z) * Radian;
  } else if (z < 0.0) {
    double const r = Sqrt(h);
    double const ri = 1.0 / r;
    return ArcTan(r * t) / r;
  } else if (abs_z < 4.7138547e-02) {
    return t * FukushimaTMaclaurin10::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 6.1227405e-02) {
    return t * FukushimaTMaclaurin11::polynomial.Evaluate(z) * Radian;
  } else if (abs_z < 7.6353468e-02) {
    return t * FukushimaTMaclaurin12::polynomial.Evaluate(z) * Radian;
  } else {
    double const r = Sqrt(-h);
    return ArcTanh(r * t) / r;
  }
}

// TODO(phl): This is extremely imprecise near large multiples of π.  Use a
// better algorithm (Payne-Hanek?).
void Reduce(Angle const& angle,
            Angle& fractional_part,
            std::int64_t& integer_part) {
  double const angle_in_half_cycles = angle / (π * Radian);
  double reduced_in_half_cycles;
#if PRINCIPIA_USE_SSE3_INTRINSICS
  auto const& x = angle_in_half_cycles;
  __m128d const x_128d = _mm_set_sd(x);
  integer_part = _mm_cvtsd_si64(x_128d);
  reduced_in_half_cycles = _mm_cvtsd_f64(
      _mm_sub_sd(x_128d,
                 _mm_cvtsi64_sd(__m128d{}, integer_part)));
#else
  integer_part = std::nearbyint(angle_in_half_cycles);
  reduced_in_half_cycles = angle_in_half_cycles - integer_part;
#endif
  fractional_part = reduced_in_half_cycles * π * Radian;
}

template<typename ThirdKind, typename>
void FukushimaEllipticBDJReduced(Angle const& φ,
                                 double const n,
                                 double const mc,
                                 Angle& B_φǀm,
                                 Angle& D_φǀm,
                                 ThirdKind& J_φ_nǀm) {
  DCHECK_LE(φ, π/2 * Radian);
  DCHECK_GE(φ, 0 * Radian);
  DCHECK_LE(0, n);
  DCHECK_GE(1, n);
  DCHECK_LE(0, mc);
  DCHECK_GE(1, mc);

  // NOTE(phl): The original Fortran code [Fuk18] had φs = 1.345 * Radian,
  // which, according to the above-mentioned paper, is suitable for single
  // precision. However, this is double precision.  Importantly, this doesn't
  // match the value of ys.  The discrepancy has a 5-10% impact on performance.
  // I am not sure if it has an impact on correctness.

  // Sin(φs)^2 must be approximately ys.
  constexpr Angle φs = 1.249 * Radian;
  constexpr double ys = 0.9;

  Angle B_m{uninitialized};        // B(m).
  Angle D_m{uninitialized};        // D(m).
  ThirdKind J_nǀm{uninitialized};  // J(n|m).

  // The selection rule in [Fuk11b] section 2.1, equations (7–11) and [Fuk12]
  // section 3.2, equations (22) and (23).  The identifiers follow Fukushima's
  // notation.
  // NOTE(phl): The computation of 1 - c² loses accuracy with respect to the
  // evaluation of Sin(φ).
  if (φ < φs) {
    FukushimaEllipticBsDsJs(Sin(φ), n, mc, B_φǀm, D_φǀm, J_φ_nǀm);
  } else {
    double const m = 1.0 - mc;
    double const nc = 1.0 - n;
    double const h = n * nc * (n - m);
    double const c = Cos(φ);
    double const c² = c * c;
    double const z²_denominator = mc + m * c²;
    if (c² < ys * z²_denominator) {
      double const z = c / Sqrt(z²_denominator);
      Angle Bs{uninitialized};      // Bs(z|m).
      Angle Ds{uninitialized};      // Ds(z|m).
      ThirdKind Js{uninitialized};  // Js(z, n|m).
      FukushimaEllipticBsDsJs(z, n, mc, Bs, Ds, Js);
      FukushimaEllipticBDJ(nc, mc, B_m, D_m, J_nǀm);
      double const sz = z * Sqrt(1.0 - c²);
      B_φǀm = B_m - (Bs - sz * Radian);
      D_φǀm = D_m - (Ds + sz * Radian);
      if constexpr (should_compute<ThirdKind>) {
        double const t = sz / nc;
        J_φ_nǀm = J_nǀm - (Js + FukushimaT(t, h));
      }
    } else {
      double const w²_numerator = mc * (1.0 - c²);
      if (w²_numerator < c² * z²_denominator) {
        FukushimaEllipticBcDcJc(c, n, mc, B_φǀm, D_φǀm, J_φ_nǀm);
      } else {
        double const w²_denominator = z²_denominator;
        double const w²_over_mc = (1.0 - c²) / w²_denominator;
        Angle Bc{uninitialized};      // Bc(w|m).
        Angle Dc{uninitialized};      // Dc(w|m).
        ThirdKind Jc{uninitialized};  // Jc(w, n|m).
        FukushimaEllipticBcDcJc(Sqrt(mc * w²_over_mc), n, mc, Bc, Dc, Jc);
        FukushimaEllipticBDJ(nc, mc, B_m, D_m, J_nǀm);
        double const sz = c * Sqrt(w²_over_mc);
        B_φǀm = B_m - (Bc - sz * Radian);
        D_φǀm = D_m - (Dc + sz * Radian);
        if constexpr (should_compute<ThirdKind>) {
          double const t = sz / nc;
          J_φ_nǀm = J_nǀm - (Jc + FukushimaT(t, h));
        }
      }
    }
  }
}

template<typename ThirdKind, typename>
void FukushimaEllipticBDJ(Angle const& φ,
                          double const n,
                          double const mc,
                          Angle& B_φǀm,
                          Angle& D_φǀm,
                          ThirdKind& J_φ_nǀm) {
  // See Appendix B of [Fuk11b] and Appendix A.1 of [Fuk12] for argument
  // reduction.

  // [Fuk12] A.1: Reduction of amplitude.
  if (φ < 0 * Radian || φ > π/2 * Radian ) {
    // TODO(phl): This is extremely imprecise near large multiples of π.  Use a
    // better algorithm (Payne-Hanek?).
    Angle φ_reduced{uninitialized};
    std::int64_t j;
    Reduce(φ, φ_reduced, j);
    Angle const abs_φ_reduced = Abs(φ_reduced);

    bool has_computed_complete_integrals = false;
    Angle B_m{uninitialized};        // B(m).
    Angle D_m{uninitialized};        // D(m).
    ThirdKind J_nǀm{uninitialized};  // J(n, m).
    FukushimaEllipticBDJ(abs_φ_reduced, n, mc, B_φǀm, D_φǀm, J_φ_nǀm);

    if (φ_reduced < 0.0 * Radian) {
      // TODO(egg): Much ado about nothing's sign bit.
      B_φǀm = -B_φǀm;
      D_φǀm = -D_φǀm;
      if constexpr (should_compute<ThirdKind>) {
        J_φ_nǀm = -J_φ_nǀm;
      }
    }
    if (j != 0) {
      double const nc = 1.0 - n;
      FukushimaEllipticBDJ(nc, mc, B_m, D_m, J_nǀm);

      // See [Fuk11b], equations (B.2), and [Fuk12], equation (A.2).
      B_φǀm += 2 * j * B_m;
      D_φǀm += 2 * j * D_m;
      if constexpr (should_compute<ThirdKind>) {
        J_φ_nǀm += 2 * j * J_nǀm;
      }
    }
    return;
  }

  // [Fuk11b] B.2, [Fuk12] A.2: Reduction of parameter.
  // NOTE(phl): Not implementing the special values mc = 0, etc. on the
  // assumption that the normal implementation will work.
  if (mc < 0) {  // m > 1
    double const m = 1 - mc;
    Angle const φR = ArcSin(Sqrt(m) * Sin(φ));
    double const nR = n / m;
    double const mR = 1 / m;
    double const mcR = -mc / m;
    double const sqrt_mR = Sqrt(mR);

    Angle B_φRǀmR{uninitialized};
    Angle D_φRǀmR{uninitialized};
    ThirdKind J_φR_nRǀmR{uninitialized};
    FukushimaEllipticBDJ(φR, nR, mcR, B_φRǀmR, D_φRǀmR, J_φR_nRǀmR);

    B_φǀm = sqrt_mR * (B_φRǀmR + mcR * D_φRǀmR);
    D_φǀm = mR * sqrt_mR * D_φRǀmR;
    if constexpr (should_compute<ThirdKind>) {
      J_φ_nǀm = mR * sqrt_mR * J_φR_nRǀmR;
    }
    return;
  } else if (mc > 1) {  // m < 0
    double const m = 1 - mc;
    double const sin_φ = Sin(φ);
    Angle const φN = ArcSin(Sqrt(mc / (1 - m * Pow<2>(sin_φ))) * sin_φ);
    double const nN = (n - m) / mc;
    double const mN = -m / mc;
    double const mcN = 1 / mc;

    Angle B_φNǀmN{uninitialized};
    Angle D_φNǀmN{uninitialized};
    ThirdKind J_φN_nNǀmN{uninitialized};
    FukushimaEllipticBDJ(φN, nN, mcN, B_φNǀmN, D_φNǀmN, J_φN_nNǀmN);

    double const sin_φN = Sin(φN);
    double const sqrt_mcN = Sqrt(mcN);
    Angle const addend =
        sin_φN * Cos(φN) / Sqrt(1 - mN * Pow<2>(sin_φN)) * Radian;
    B_φǀm = sqrt_mcN * (D_φNǀmN + addend);
    D_φǀm = sqrt_mcN * (B_φNǀmN - addend);
    if constexpr (should_compute<ThirdKind>) {
      J_φ_nǀm = mN * Sqrt(mN) * J_φN_nNǀmN;
    }
    return;
  }

  // [Fuk12] A.3: Reduction of characteristics.
  // NOTE(phl): Same as above, not implementing the special values.
  if constexpr (should_compute<ThirdKind>) {
    if (n > 1) {
      double const m = 1 - mc;
      double const nc = 1 - n;
      double const t1 = Tan(φ) / Sqrt(1 - m * Pow<2>(Sin(φ)));
      double const h1 = nc * (n - m) / n;
      double const n1 = m / n;

      ThirdKind J_φ_n1ǀm{uninitialized};
      FukushimaEllipticBDJ(φ, n1, mc, B_φǀm, D_φǀm, J_φ_n1ǀm);

      J_φ_nǀm = (-B_φǀm - D_φǀm + FukushimaT(t1, h1) - n1 * J_φ_n1ǀm) / n;
      return;
    } else if (n < 0) {
      double const m = 1 - mc;
      double const nc = 1 - n;
      double const sin_φ = Sin(φ);
      double const t2 = sin_φ * Cos(φ) / Sqrt(1 - m * Pow<2>(sin_φ));
      double const h2 = -n * (m - n) / (1 - n);
      double const n2 = (m - n) / nc;

      ThirdKind J_φ_n2ǀm{uninitialized};
      FukushimaEllipticBDJ(φ, n2, mc, B_φǀm, D_φǀm, J_φ_n2ǀm);

      J_φ_nǀm = (B_φǀm + D_φǀm + FukushimaT(t2, h2) - n2 * J_φ_n2ǀm) / nc;
      return;
    }
  }

  // No further reduction needed.
  FukushimaEllipticBDJReduced(φ, n, mc, B_φǀm, D_φǀm, J_φ_nǀm);
}

template<typename ThirdKind, typename>
void EllipticFEΠ(Angle const& φ,
                 double const n,
                 double const mc,
                 Angle& F_φǀm,
                 Angle& E_φǀm,
                 ThirdKind& Π_φ_nǀm) {
  Angle B{uninitialized};
  Angle D{uninitialized};
  ThirdKind J{uninitialized};
  FukushimaEllipticBDJ(φ, n, mc, B, D, J);
  F_φǀm = B + D;
  E_φǀm = B + mc * D;
  if constexpr (should_compute<ThirdKind>) {
    Π_φ_nǀm = F_φǀm + n * J;
  }
}

}  // namespace

void FukushimaEllipticBDJ(Angle const& φ,
                          double const n,
                          double const mc,
                          Angle& B_φǀm,
                          Angle& D_φǀm,
                          Angle& J_φ_nǀm) {
  return FukushimaEllipticBDJ<Angle>(φ, n, mc, B_φǀm, D_φǀm, J_φ_nǀm);
}

void FukushimaEllipticBD(Angle const& φ,
                         double const mc,
                         Angle& B_φǀm,
                         Angle& D_φǀm) {
  FukushimaEllipticBDJ(φ, /*n=*/1, mc, B_φǀm, D_φǀm, /*J_φ_nǀm=*/unused);
}

Angle EllipticF(Angle const& φ, double const mc) {
  Angle F{uninitialized};
  Angle E{uninitialized};
  EllipticFEΠ(φ, /*n=*/1, mc, F, E, /*Π=*/unused);
  return F;
}

Angle EllipticE(Angle const& φ, double const mc) {
  Angle F{uninitialized};
  Angle E{uninitialized};
  EllipticFEΠ(φ, /*n=*/1, mc, F, E, /*Π=*/unused);
  return E;
}

Angle EllipticΠ(Angle const& φ, double const n, double const mc) {
  Angle F{uninitialized};
  Angle E{uninitialized};
  Angle Π{uninitialized};
  EllipticFEΠ(φ, n, mc, F, E, Π);
  return Π;
}

void EllipticFE(Angle const& φ, double mc, Angle& F_φǀm, Angle& E_φǀm) {
  EllipticFEΠ(φ, /*n=*/1, mc, F_φǀm, E_φǀm, /*Π=*/unused);
}

void EllipticFEΠ(Angle const& φ,
                 double const n,
                 double const mc,
                 Angle& F_φǀm,
                 Angle& E_φǀm,
                 Angle& Π_φ_nǀm) {
  EllipticFEΠ<Angle>(φ, n, mc, F_φǀm, E_φǀm, Π_φ_nǀm);
}

// Note that the identifiers in the function definition are not the same as
// those in the function declaration.
// The notation here follows [Fuk09], whereas the notation in the function
// declaration uses |mc| for consistency with the other functions.
Angle EllipticK(double const mʹ) {
  DCHECK_LE(0, mʹ);
  DCHECK_GE(1, mʹ);
  // TODO(phl): Use a binary split of [0, 1] to reduce the number of
  // comparisons.
  double const m = 1.0 - mʹ;
  if (m == 0.0) {
    return π / 2 * Radian;
  } else if (mʹ < std::numeric_limits<double>::epsilon() / 2.0) {
    return (2.0 * log_2 - 0.5 * std::log(mʹ)) * Radian;
  } else if (mʹ < 0.1) {
    // The complementary nome.
    double const qʹ = EllipticNomeQ<14>(mʹ);
    // Use K′ = K(m′), see [Fuk09], equations (15) and (29).
    double const Kʹ = elliptic_K_taylor_0_05.Evaluate(mʹ - 0.05);
    return -Kʹ * (1 / π) * std::log(qʹ) * Radian;
  } else if (m <= 0.1) {
    return elliptic_K_taylor_0_05.Evaluate(m - 0.05) * Radian;
  } else if (m <= 0.2) {
    return elliptic_K_taylor_0_15.Evaluate(m - 0.15) * Radian;
  } else if (m <= 0.3) {
    return elliptic_K_taylor_0_25.Evaluate(m - 0.25) * Radian;
  } else if (m <= 0.4) {
    return elliptic_K_taylor_0_35.Evaluate(m - 0.35) * Radian;
  } else if (m <= 0.5) {
    return elliptic_K_taylor_0_45.Evaluate(m - 0.45) * Radian;
  } else if (m <= 0.6) {
    return elliptic_K_taylor_0_55.Evaluate(m - 0.55) * Radian;
  } else if (m <= 0.7) {
    return elliptic_K_taylor_0_65.Evaluate(m - 0.65) * Radian;
  } else if (m <= 0.8) {
    return elliptic_K_taylor_0_75.Evaluate(m - 0.75) * Radian;
  } else if (m <= 0.85) {
    return elliptic_K_taylor_0_825.Evaluate(m - 0.825) * Radian;
  } else {
    return elliptic_K_taylor_0_875.Evaluate(m - 0.875) * Radian;
  }
}

}  // namespace internal_elliptic_integrals
}  // namespace numerics
}  // namespace principia
