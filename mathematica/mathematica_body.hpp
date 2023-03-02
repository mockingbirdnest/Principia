#pragma once

#include "mathematica/mathematica.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "base/mod.hpp"
#include "base/not_constructible.hpp"
#include "base/traits.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace mathematica {
namespace internal_mathematica {

using astronomy::J2000;
using namespace principia::base::_mod;
using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::quantities::_quantities;

// Wraps the string in quotes and escapes things properly.
inline std::string Escape(std::string_view const str) {
  std::string result = "\"";
  for (const char c : str) {
    switch (c) {
      case '"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
        break;
      default:
        result += c;
        break;
    }
  }
  result += "\"";
  return result;
}

// Does not wrap its arguments in ToMathematica.
inline std::string RawApply(
    std::string const& function,
    std::vector<std::string> const& arguments) {
  std::string result;
  result += function;
  result += "[";
  for (int i = 0; i < arguments.size(); ++i) {
    result += arguments[i];
    result += (i + 1 < arguments.size() ? "," : "");
  }
  result += "]";
  return result;
}

// A helper struct to scan the elements of a tuple and stringify them.
template<int index, typename Tuple, typename OptionalExpressIn>
struct TupleHelper : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   OptionalExpressIn express_in) {
    TupleHelper<index - 1, Tuple, OptionalExpressIn>::ToMathematicaStrings(
        tuple, expressions, express_in);
    expressions.push_back(ToMathematica(std::get<index - 1>(tuple),
                          express_in));
  }
};

template<typename Tuple, typename OptionalExpressIn>
struct TupleHelper<0, Tuple, OptionalExpressIn> : not_constructible {
  static void ToMathematicaStrings(Tuple const& tuple,
                                   std::vector<std::string>& expressions,
                                   OptionalExpressIn express_in) {}
};

template<typename V, typename A, int d,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematicaBody(
    PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
    OptionalExpressIn express_in) {
  using Coefficients =
      typename PolynomialInMonomialBasis<V, A, d, E>::Coefficients;
  std::vector<std::string> coefficients;
  coefficients.reserve(std::tuple_size_v<Coefficients>);
  TupleHelper<std::tuple_size_v<Coefficients>,
              Coefficients,
              OptionalExpressIn>::ToMathematicaStrings(polynomial.coefficients_,
                                                       coefficients,
                                                       express_in);
  std::string argument;
  if constexpr (is_instance_of_v<Point, A>) {
    argument = RawApply("Subtract",
                        {"#", ToMathematica(polynomial.origin_, express_in)});
  } else {
    argument = "#";
  }
  std::vector<std::string> monomials;
  for (int i = 0; i < coefficients.size(); ++i) {
    if (i == 0) {
      monomials.push_back(coefficients[i]);
    } else if (i == 1) {
      monomials.push_back(RawApply("Times", {coefficients[i], argument}));
    } else {
      monomials.push_back(RawApply(
          "Times",
          {coefficients[i], RawApply("Power", {argument, std::to_string(i)})}));
    }
  }
  return RawApply("Plus", monomials);
}

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematicaBody(
    PoissonSeries<V, ad, pd, E> const& series,
    OptionalExpressIn express_in) {
  std::vector<std::string> components = {
      ToMathematicaBody(series.aperiodic_, express_in)};
  for (auto const& [ω, polynomials] : series.periodic_) {
    std::string const polynomial_sin =
        ToMathematicaBody(polynomials.sin, express_in);
    std::string const polynomial_cos =
        ToMathematicaBody(polynomials.cos, express_in);
    std::string const angle =
        RawApply("Times",
                 {ToMathematica(ω, express_in),
                  RawApply("Subtract",
                           {"#", ToMathematica(series.origin_, express_in)})});
    components.push_back(
        RawApply("Times", {polynomial_sin, RawApply("Sin", {angle})}));
    components.push_back(
        RawApply("Times", {polynomial_cos, RawApply("Cos", {angle})}));
  }
  return RawApply("Plus", components);
}

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematicaBody(
    PiecewisePoissonSeries<V, ad, pd, E> const& series,
    OptionalExpressIn express_in) {
  std::vector<std::string> conditions_and_functions;
  for (int i = 0; i < series.series_.size(); ++i) {
    std::string const function =
        ToMathematicaBody(series.series_[i], express_in);
    std::string const condition = RawApply(
        "Between",
        {"#",
         RawApply("List",
                  {ToMathematica(series.bounds_[i], express_in),
                   ToMathematica(series.bounds_[i + 1], express_in)})});
    conditions_and_functions.push_back(RawApply("List", {function, condition}));
  }
  return RawApply("Piecewise", {RawApply("List", conditions_and_functions)});
}

template<typename... Qs>
template<typename Q>
double ExpressIn<Qs...>::operator()(Q const& q) const {
  return Divide<Q::Dimensions::Length, Length>(
      Divide<Q::Dimensions::Mass, Mass>(
          Divide<Q::Dimensions::Time, Time>(
              Divide<Q::Dimensions::Current, Current>(
                  Divide<Q::Dimensions::Temperature, Temperature>(
                      Divide<Q::Dimensions::Amount, Amount>(
                          Divide<Q::Dimensions::LuminousIntensity,
                                 LuminousIntensity>(
                              Divide<Q::Dimensions::Angle, Angle>(q))))))));
}

template<typename... Qs>
template<std::int64_t exponent, typename Q1, typename Q2>
Quotient<Q2, Exponentiation<Q1, exponent>> ExpressIn<Qs...>::Divide(
    Q2 const& q2) const {
  if constexpr (exponent == 0) {
    return q2;
  } else {
    return q2 / Pow<exponent>(std::get<Q1>(units_));
  }
}

template<typename T, typename OptionalExpressIn>
std::string Apply(std::string const& name,
                  T const& right,
                  OptionalExpressIn express_in) {
  return RawApply("Apply", {name, ToMathematica(right, express_in)});
}

inline std::string Evaluate(std::string const& expression) {
  return RawApply("Evaluate", {expression});
}

template<typename T, typename OptionalExpressIn>
std::string Rule(std::string const& name,
                 T const& right,
                 OptionalExpressIn express_in) {
  return RawApply("Rule", {name, ToMathematica(right, express_in)});
}

template<typename T, typename OptionalExpressIn>
std::string Set(std::string const& name,
                T const& right,
                OptionalExpressIn express_in) {
  return RawApply("Set", {name, ToMathematica(right, express_in)}) + ";\n";
}

template<typename T, typename U, typename OptionalExpressIn>
std::string PlottableDataset(std::vector<T> const& x,
                             std::vector<U> const& y,
                             OptionalExpressIn express_in) {
  std::vector<std::string> const xy = {ToMathematica(x, express_in),
                                       ToMathematica(y, express_in)};
  return RawApply("Transpose", {RawApply("List", xy)});
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(std::vector<T> const& list,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(list.size());
  for (auto const& expression : list) {
    expressions.emplace_back(ToMathematica(expression, express_in));
  }
  return RawApply("List", expressions);
}

template<typename It, typename OptionalExpressIn>
std::string ToMathematica(It const begin, It const end,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  for (auto it = begin; it != end; ++it) {
    expressions.emplace_back(ToMathematica(*it, express_in));
  }
  return RawApply("List", expressions);
}

template<typename OptionalExpressIn>
std::string ToMathematica(bool const b, OptionalExpressIn /*express_in*/) {
  return b ? "True" : "False";
}

template<typename T, typename, typename OptionalExpressIn>
std::string ToMathematica(T const integer, OptionalExpressIn /*express_in*/) {
  return std::to_string(integer);
}

template<typename T, typename, typename OptionalExpressIn, typename>
std::string ToMathematica(T const real,
                          OptionalExpressIn /*express_in*/) {
  std::string absolute_value;
  if (std::isinf(real)) {
    absolute_value = "Infinity";
  } else if (std::isnan(real)) {
    absolute_value = "Indeterminate";
  } else if (real == 0) {
    absolute_value = "0";
  } else {
    constexpr int τ = std::numeric_limits<T>::digits;
    int const exponent = std::ilogb(real);
    // This offset makes n an integer in [β^(τ-1), β^τ[, i.e., a τ-digit
    // integer.
    int exponent_offset = τ - 1;
    if (std::numeric_limits<T>::radix == 2) {
      // For binary floating point, push the leading 1 to the least significant
      // bit of a hex digit.
      exponent_offset += mod(1 - τ, 4);
    }
    std::int64_t const n =
        std::scalbln(std::abs(real), exponent_offset - exponent);
    absolute_value =
        RawApply("Times",
                 {std::numeric_limits<T>::radix == 10
                      ? ToMathematica(n)
                      : (std::stringstream()
                         << "16^^" << std::uppercase << std::hex << n)
                            .str(),
                  RawApply("Power",
                           {"2",
                            RawApply("Subtract",
                                     {ToMathematica(std::ilogb(real)),
                                      ToMathematica(exponent_offset)})})});
  }
  return std::signbit(real) ? RawApply("Minus", {absolute_value})
                            : absolute_value;
}

template<typename OptionalExpressIn>
std::string ToMathematica(Quaternion const& quaternion,
                          OptionalExpressIn /*express_in*/) {
  return RawApply("Quaternion",
                  {ToMathematica(quaternion.real_part()),
                   ToMathematica(quaternion.imaginary_part().x),
                   ToMathematica(quaternion.imaginary_part().y),
                   ToMathematica(quaternion.imaginary_part().z)});
}

template<typename T, int size, typename OptionalExpressIn>
std::string ToMathematica(FixedVector<T, size> const& fixed_vector,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  for (int i = 0; i < size; ++i) {
    expressions.emplace_back(ToMathematica(fixed_vector[i], express_in));
  }
  return RawApply("List", expressions);
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(R3Element<T> const& r3_element,
                          OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.emplace_back(ToMathematica(r3_element.x, express_in));
  expressions.emplace_back(ToMathematica(r3_element.y, express_in));
  expressions.emplace_back(ToMathematica(r3_element.z, express_in));
  return RawApply("List", expressions);
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(R3x3Matrix<T> const& r3x3_matrix,
                          OptionalExpressIn express_in) {
  std::vector<R3Element<T>> rows;
  rows.push_back(r3x3_matrix.row_x());
  rows.push_back(r3x3_matrix.row_y());
  rows.push_back(r3x3_matrix.row_z());
  return ToMathematica(rows, express_in);
}

template<typename D, typename... Qs>
std::string ToMathematica(Quantity<D> const& quantity,
                          ExpressIn<Qs...> express_in) {
  return ToMathematica(express_in(quantity));
}

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity,
                          decltype(PreserveUnits) express_in) {
  std::string s = DebugString(quantity);
  std::string const number = ToMathematica(quantity / si::Unit<Quantity<D>>);
  std::size_t const split = s.find(" ");
  std::string const units = Escape(s.substr(split, s.size()));
  return RawApply("Quantity", {number, units});
}

template<typename D>
std::string ToMathematica(Quantity<D> const& quantity,
                          std::nullopt_t express_in) {
  static_assert(
      !std::is_same_v<Quantity<D>, Quantity<D>>,
      "Must specify a way to express units for dimensionful quantities");
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(DoublePrecision<T> const& double_precision,
                          OptionalExpressIn express_in) {
  return RawApply("Plus",
                  {ToMathematica(double_precision.value, express_in),
                   ToMathematica(double_precision.error, express_in)});
}

template<typename S, typename F, typename OptionalExpressIn>
std::string ToMathematica(Vector<S, F> const& vector,
                          OptionalExpressIn express_in) {
  return ToMathematica(vector.coordinates(), express_in);
}

template<typename S, typename F, typename OptionalExpressIn>
std::string ToMathematica(Bivector<S, F> const& bivector,
                          OptionalExpressIn express_in) {
  return ToMathematica(bivector.coordinates(), express_in);
}

template<typename V, typename OptionalExpressIn>
std::string ToMathematica(Point<V> const & point,
                          OptionalExpressIn express_in) {
  return ToMathematica(point - Point<V>(), express_in);
}

template<typename S, typename F,
         template<typename, typename> typename M,
         typename OptionalExpressIn>
std::string ToMathematica(SymmetricBilinearForm<S, F, M> const& form,
                          OptionalExpressIn express_in) {
  return ToMathematica(form.coordinates(), express_in);
}

template<typename F, typename OptionalExpressIn>
std::string ToMathematica(DegreesOfFreedom<F> const& degrees_of_freedom,
                          OptionalExpressIn express_in) {
  return RawApply(
      "List",
      std::vector<std::string>{
          ToMathematica(degrees_of_freedom.position(), express_in),
          ToMathematica(degrees_of_freedom.velocity(), express_in)});
}

template<typename Tuple, typename, typename OptionalExpressIn>
std::string ToMathematica(Tuple const& tuple, OptionalExpressIn express_in) {
  std::vector<std::string> expressions;
  expressions.reserve(std::tuple_size_v<Tuple>);
  TupleHelper<std::tuple_size_v<Tuple>, Tuple, OptionalExpressIn>::
      ToMathematicaStrings(tuple, expressions, express_in);
  return RawApply("List", expressions);
}

template<typename Scalar, typename OptionalExpressIn>
std::string ToMathematica(UnboundedLowerTriangularMatrix<Scalar> const& matrix,
                          OptionalExpressIn express_in) {
  std::vector<std::string> rows;
  rows.reserve(matrix.rows());
  for (int i = 0; i < matrix.rows(); ++i) {
    std::vector<std::string> row;
    row.reserve(matrix.rows());
    for (int j = 0; j <= i; ++j) {
      row.push_back(ToMathematica(matrix(i, j), express_in));
    }
    for (int j = i + 1; j < matrix.rows(); ++j) {
      row.push_back(ToMathematica(Scalar{}, express_in));
    }
    rows.push_back(RawApply("List", row));
  }
  return RawApply("List", rows);
}

template<typename Scalar, typename OptionalExpressIn>
std::string ToMathematica(UnboundedUpperTriangularMatrix<Scalar> const& matrix,
                          OptionalExpressIn express_in) {
  std::vector<std::string> rows;
  rows.reserve(matrix.columns());
  for (int i = 0; i < matrix.columns(); ++i) {
    std::vector<std::string> row;
    row.reserve(matrix.columns());
    for (int j = 0; j < i; ++j) {
      row.push_back(ToMathematica(Scalar{}, express_in));
    }
    for (int j = i; j < matrix.columns(); ++j) {
      row.push_back(ToMathematica(matrix(i, j), express_in));
    }
    rows.push_back(RawApply("List", row));
  }
  return RawApply("List", rows);
}

template<typename Scalar, typename OptionalExpressIn>
std::string ToMathematica(UnboundedVector<Scalar> const& vector,
                          OptionalExpressIn express_in) {
  std::vector<std::string> elements;
  elements.reserve(vector.size());
  for (int i = 0; i < vector.size(); ++i) {
    elements.push_back(ToMathematica(vector[i], express_in));
  }
  return RawApply("List", elements);
}

template<typename R, typename, typename, typename OptionalExpressIn>
std::string ToMathematica(R const ref,
                          OptionalExpressIn express_in) {
  return RawApply("List",
                  std::vector<std::string>{
                      ToMathematica(ref.time, express_in),
                      ToMathematica(ref.degrees_of_freedom, express_in)});
}

template<typename V, typename A, int d,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematica(
    PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
    OptionalExpressIn express_in) {
  return RawApply("Function", {ToMathematicaBody(polynomial, express_in)});
}

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematica(PoissonSeries<V, ad, pd, E> const& series,
                          OptionalExpressIn express_in) {
  return RawApply("Function", {ToMathematicaBody(series, express_in)});
}

template<typename V, int ad, int pd,
         template<typename, typename, int> class E,
         typename OptionalExpressIn>
std::string ToMathematica(PiecewisePoissonSeries<V, ad, pd, E> const& series,
                          OptionalExpressIn express_in) {
  return RawApply("Function", {ToMathematicaBody(series, express_in)});
}

template<typename OptionalExpressIn>
std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements,
    OptionalExpressIn express_in) {
  return ToMathematica(std::make_tuple((elements.t - J2000),
                                       elements.a,
                                       elements.h,
                                       elements.k,
                                       elements.λ,
                                       elements.p,
                                       elements.q,
                                       elements.pʹ,
                                       elements.qʹ),
                       express_in);
}

template<typename T, typename OptionalExpressIn>
std::string ToMathematica(std::optional<T> const& opt,
                          OptionalExpressIn express_in) {
  std::vector<T> value;
  if (opt.has_value()) {
    value.push_back(opt.value());
  }
  return ToMathematica(value, express_in);
}

template<typename OptionalExpressIn>
std::string ToMathematica(char const* const str,
                          OptionalExpressIn /*express_in*/) {
  return Escape(str);
}

template<typename OptionalExpressIn>
std::string ToMathematica(std::string const& str,
                          OptionalExpressIn /*express_in*/) {
  return Escape(str);
}

}  // namespace internal_mathematica
}  // namespace mathematica
}  // namespace principia
