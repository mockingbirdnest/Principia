#pragma once

#include <optional>
#include <string>

#include "absl/strings/str_replace.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

using quantities::Abs;
using quantities::Format;
using quantities::SIUnit;

namespace si = quantities::si;

template<typename Dimensions>
Quantity<Dimensions> ApproximateQuantity<Quantity<Dimensions>>::min() const {
  return min_multiplier_ * unit_;
}

template<typename Dimensions>
Quantity<Dimensions> ApproximateQuantity<Quantity<Dimensions>>::max() const {
  return max_multiplier_ * unit_;
}

template<typename Dimensions>
Quantity<Dimensions> ApproximateQuantity<Quantity<Dimensions>>::unit() const {
  return unit_;
}

template<typename Dimensions>
bool ApproximateQuantity<Quantity<Dimensions>>::has_trivial_unit() const {
  return unit_ == SIUnit<Quantity<Dimensions>>();
}

template<typename Dimensions>
double ApproximateQuantity<Quantity<Dimensions>>::UlpDistance(
    Quantity<Dimensions> const& q) const {
  return ulp_ * Abs(q - (min_multiplier_ + max_multiplier_) * 0.5 * unit_) /
         ((max_multiplier_ - min_multiplier_) * 0.5 * unit_);
}

template<typename Dimensions>
std::string ApproximateQuantity<Quantity<Dimensions>>::DebugString() const {
  if (has_trivial_unit()) {
    return (negated_ ? "-" : "") + representation_ + "(" +
           std::to_string(ulp_) + ") " + Format<Dimensions>();
  } else {
    return (negated_ ? "-" : "") + representation_ +
           "(" + std::to_string(ulp_) + ") * " + quantities::DebugString(unit_);
  }
}

template<typename Dimensions>
ApproximateQuantity<Quantity<Dimensions>>::ApproximateQuantity(
    std::string const& representation,
    int const ulp,
    bool const negated,
    double const min_multiplier,
    double const max_multiplier,
    Quantity<Dimensions> const& unit)
    : representation_(representation),
      ulp_(ulp),
      negated_(negated),
      min_multiplier_(min_multiplier),
      max_multiplier_(max_multiplier),
      unit_(unit) {}

inline ApproximateQuantity<double> ApproximateQuantity<double>::Parse(
    std::string_view const representation,
    int const ulp) {
  std::string error_representation(representation);
  std::optional<int> last_digit_index;
  bool const is_hexadecimal =
      error_representation.size() >= 2 &&
      error_representation[0] == '0' &&
      (error_representation[1] == 'x' || error_representation[1] == 'X');

  // Replace all the digits before the exponent by zeroes, except for the last
  // one which get the number of ulps.  The result is the string representation
  // of the error on the quantity.
  for (int i = 0; i < error_representation.size(); ++i) {
    char const c = error_representation[i];
    if (c >= '0' && c <= '9') {
      error_representation[i] = '0';
      last_digit_index = i;
    } else if (is_hexadecimal &&
               ((c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'))) {
      error_representation[i] = '0';
      last_digit_index = i;
    } else if ((!is_hexadecimal && (c == 'e' || c == 'E')) ||
               (is_hexadecimal && (c == 'p' || c == 'P'))) {
      CHECK(last_digit_index.has_value());
      break;
    }
  }
  if (ulp <= 9) {
    error_representation[*last_digit_index] = '0' + ulp;
  } else {
    CHECK(is_hexadecimal);
    error_representation[*last_digit_index] = 'A' + ulp - 10;
  }

  // Apparently the stupid language doesn't know how to parse literals with
  // quotes...
  std::string const stripped_representation =
      absl::StrReplaceAll(representation, {{"'", ""}});
  std::string const stripped_error_representation =
      absl::StrReplaceAll(error_representation, {{"'", ""}});
  double const value = std::strtod(stripped_representation.data(),
                                   nullptr);
  double const error = std::strtod(stripped_error_representation.c_str(),
                                   nullptr);
  return ApproximateQuantity<double>(representation,
                                     ulp,
                                     /*negated=*/false,
                                     value - error,
                                     value + error);
}

inline double ApproximateQuantity<double>::min() const {
  return min_multiplier_;
}

inline double ApproximateQuantity<double>::max() const {
  return max_multiplier_;
}

inline double ApproximateQuantity<double>::UlpDistance(double const d) const {
  return ulp_ * Abs(d - (min_multiplier_ + max_multiplier_) * 0.5) /
         ((max_multiplier_ - min_multiplier_) * 0.5);
}

inline std::string ApproximateQuantity<double>::DebugString() const {
  return (negated_ ? "-" : "") + representation_ +
         "(" + std::to_string(ulp_) + ")";
}

inline ApproximateQuantity<double>::ApproximateQuantity(
    std::string_view const representation,
    int const ulp,
    bool const negated,
    double const min_multiplier,
    double const max_multiplier)
    : representation_(representation),
      ulp_(ulp),
      negated_(negated),
      min_multiplier_(min_multiplier),
      max_multiplier_(max_multiplier) {}

template<typename Right>
ApproximateQuantity<Right> operator+(ApproximateQuantity<Right> const& right) {
  return right;
}

template<typename Right>
ApproximateQuantity<Right> operator-(ApproximateQuantity<Right> const& right) {
  return ApproximateQuantity<Right>(right.representation_,
                                    right.ulp_,
                                    !right.negated_,
                                    -right.max_multiplier_,
                                    -right.min_multiplier_,
                                    right.unit_);
}

template<>
inline ApproximateQuantity<double> operator-(
    ApproximateQuantity<double> const& right) {
  return ApproximateQuantity<double>(right.representation_,
                                     right.ulp_,
                                     !right.negated_,
                                     -right.max_multiplier_,
                                     -right.min_multiplier_);
}

template<typename Left, typename RDimensions>
ApproximateQuantity<Product<Left, Quantity<RDimensions>>> operator*(
    ApproximateQuantity<Left> const& left,
    Quantity<RDimensions> const& right) {
  return ApproximateQuantity<Product<Left, Quantity<RDimensions>>>(
      left.representation_,
      left.ulp_,
      left.negated_,
      left.min_multiplier_,
      left.max_multiplier_,
      left.unit_ * right);
}

template<typename Left, typename RDimensions>
ApproximateQuantity<Quotient<Left, Quantity<RDimensions>>> operator/(
    ApproximateQuantity<Left> const& left,
    Quantity<RDimensions> const& right) {
  return ApproximateQuantity<Quotient<Left, Quantity<RDimensions>>>(
      left.representation_,
      left.ulp_,
      left.negated_,
      left.min_multiplier_,
      left.max_multiplier_,
      left.unit_ / right);
}

template<typename Quantity>
std::ostream& operator<<(std::ostream& out,
                         ApproximateQuantity<Quantity> const& q) {
  out << q.DebugString();
  return out;
}

#define PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(symbol, ulp) \
  inline ApproximateQuantity<double> operator""_##symbol(               \
      char const* const representation) {                               \
    return ApproximateQuantity<double>::Parse(representation, ulp);     \
  }

PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑴, 1)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑵, 2)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑶, 3)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑷, 4)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑸, 5)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑹, 6)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑺, 7)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑻, 8)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(⑼, 9)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄐, 0xA)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄑, 0xB)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄒, 0xC)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄓, 0xD)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄔, 0xE)
PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION(🄕, 0xF)

#undef PRINCIPIA_APPROXIMATE_QUANTITY_OPERATOR_DEFINITION

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia
