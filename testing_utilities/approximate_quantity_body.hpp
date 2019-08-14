#pragma once

#include <optional>

#include "testing_utilities/approximate_quantity.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

template<typename Quantity>
Quantity ApproximateQuantity<Quantity>::min() const {
  return min_multiplier_ * unit_;
}

template<typename Quantity>
Quantity ApproximateQuantity<Quantity>::max() const {
  return max_multiplier_ * unit_;
}

template<typename Quantity>
ApproximateQuantity<Quantity>::ApproximateQuantity(
    std::string const& representation,
    Quantity const& unit,
    double const min_multiplier,
    double const max_multiplier)
    : representation_(representation),
      unit_(unit),
      min_multiplier_(min_multiplier),
      max_multiplier_(max_multiplier) {}

template<typename Left, typename Right>
Product<Left, Right> operator*(ApproximateQuantity<Left> const& left,
                               Right const& right) {
  return ApproximateQuantity<Product<Left, Right>>(left.representation_,
                                                   left.unit_ * right,
                                                   left.min_multiplier_,
                                                   left.max_multiplier_);
}

ApproximateQuantity<double> operator""_⑴(char const* representation) {
  std::string error_representation(representation);
  std::optional<int> last_digit_index;
  bool const is_hexadecimal =
      error_representation.size() >= 2 &&
      error_representation[0] == '0' &&
      (error_representation[1] == 'x' || error_representation[1] == 'X');
  for (int i = 0; i < error_representation.size(); ++i) {
    char const c = error_representation[i];
    if (c >= '1' && c <= '9') {
      error_representation[i] = '0';
      last_digit_index = i;
    } else if (is_hexadecimal &&
               ((c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'))) {
      error_representation[i] = '0';
      last_digit_index = i;
    } else if ((!is_hexadecimal && (c == 'e' || c == 'E')) ||
               (is_hexadecimal && (c == 'p' || c == 'P'))) {
      CHECK(last_digit_index);
      break;
    }
  }
  error_representation[*last_digit_index] = '1';
  double const value = std::strtod(representation, nullptr);
  double const error = std::strtod(error_representation.c_str(), nullptr);
  LOG(ERROR)<<error_representation;
  LOG(ERROR)<<error;
  return ApproximateQuantity<double>(representation,
                                     /*unit=*/1.0,
                                     value - error,
                                     value + error);
}

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia
