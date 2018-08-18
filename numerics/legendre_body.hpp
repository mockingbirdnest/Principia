
#pragma once

#include "numerics/legendre.hpp"

#include "base/not_constructible.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

using base::not_constructible;

//TODO(phl): Do I *really* have to do this?
template<typename Tuple, int k = 0, int size = std::tuple_size_v<Tuple>>
struct TupleArithmetic : not_constructible {
  static constexpr Tuple Add(Tuple const& left, Tuple const& right);
  static constexpr Tuple Multiply(int left, Tuple const& right);
};

template<typename Tuple, int size>
struct TupleArithmetic<Tuple, size, size> : not_constructible {
  static constexpr Tuple Add(Tuple const& left, Tuple const& right);
  static constexpr Tuple Multiply(int left, Tuple const& right);
};

template<typename Tuple, int k, int size>
constexpr Tuple TupleArithmetic<Tuple, k, size>::Add(Tuple const& left,
                                                     Tuple const& right) {
  Tuple result = right;
  std::get<k>(result) += std::get<k>(left);
  return TupleArithmetic<Tuple, k + 1, size>::Add(left, result);
}

template<typename Tuple, int k, int size>
constexpr Tuple TupleArithmetic<Tuple, k, size>::Multiply(int const left,
                                                          Tuple const& right) {
  Tuple result = right;
  std::get<k>(result) *= left;
  return TupleArithmetic<Tuple, k + 1, size>::Multiply(left, result);
}

template<typename Tuple, int size>
constexpr Tuple TupleArithmetic<Tuple, size, size>::Add(Tuple const& left,
                                                        Tuple const& right) {
  return right;
}

template<typename Tuple, int size>
constexpr Tuple TupleArithmetic<Tuple, size, size>::Multiply(
    int const left,
    Tuple const& right) {
  return right;
}

template<int degree_>
struct CoefficientsGenerator {
  using Type = internal_polynomial::
      Derivatives<double, double, std::make_integer_sequence<int, degree_ + 1>>;
  static constexpr auto a =
      TupleArithmetic<typename CoefficientsGenerator<degree_ - 2>::Type>::
          Multiply(degree_ - 1, CoefficientsGenerator<degree_ - 2>::value);
  static constexpr auto b = std::tuple_cat(std::make_pair(0.0, 0.0), a);
  static constexpr auto c =
      TupleArithmetic<typename CoefficientsGenerator<degree_ - 1>::Type>::
          Multiply(2 * degree_ - 1, CoefficientsGenerator<degree_ - 1>::value);
  static constexpr auto d = std::tuple_cat(c, std::make_tuple(0.0));
  static constexpr auto e = TupleArithmetic<Type>::Add(b, d);
  static constexpr Type value = e;
};

template<>
struct CoefficientsGenerator<0> {
  using Type = internal_polynomial::
      Derivatives<double, double, std::make_integer_sequence<int, 1>>;
  static constexpr Type value{1};
};

template<>
struct CoefficientsGenerator<1> {
  using Type = internal_polynomial::
      Derivatives<double, double, std::make_integer_sequence<int, 2>>;
  static constexpr Type value{0, 1};
};

template<int degree_, template<typename, typename, int> class Evaluator>
LegendrePolynomial<degree_, Evaluator>::LegendrePolynomial()
    : PolynomialInMonomialBasis(CoefficientsGenerator<degree_>::value) {}

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia
