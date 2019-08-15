#pragma once

#include <string>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

using quantities::Product;
using quantities::Quantity;

template<typename Quantity>
class ApproximateQuantity;

template<typename Dimensions>
class ApproximateQuantity<Quantity<Dimensions>> {
 public:
  Quantity<Dimensions> min() const;
  Quantity<Dimensions> max() const;

 private:
  ApproximateQuantity(std::string const& representation,
                      Quantity<Dimensions> const& unit,
                      double min_multiplier,
                      double max_multiplier);

  // The original representation.
  std::string representation_;

  // The unit for the approximate quantity.  Could for instance be Degree for an
  // angle.
  Quantity<Dimensions> unit_;

  // The interval for the approximate quantity, expressed as multiples of unit_.
  double min_multiplier_;
  double max_multiplier_;

  template<typename Left, typename Right>
  friend Product<Left, Right> operator*(ApproximateQuantity<Left> const& left,
                                        Right const& right);
};

template<>
class ApproximateQuantity<double> {
 public:
  static ApproximateQuantity<double> Parse(char const* representation, int ulp);

  double min() const;
  double max() const;

 private:
  ApproximateQuantity(std::string const& representation,
                      double min_multiplier,
                      double max_multiplier);

  // The original representation.
  std::string representation_;

  // The interval for the approximate quantity.
  double min_multiplier_;
  double max_multiplier_;
};

template<typename Left, typename Right>
Product<Left, Right> operator*(ApproximateQuantity<Left> const& left,
                               Right const& right);

ApproximateQuantity<double> operator""_⑴(char const* representation);
ApproximateQuantity<double> operator""_⑵(char const* representation);
ApproximateQuantity<double> operator""_⑶(char const* representation);
ApproximateQuantity<double> operator""_⑷(char const* representation);
ApproximateQuantity<double> operator""_⑸(char const* representation);
ApproximateQuantity<double> operator""_⑹(char const* representation);
ApproximateQuantity<double> operator""_⑺(char const* representation);
ApproximateQuantity<double> operator""_⑻(char const* representation);
ApproximateQuantity<double> operator""_⑼(char const* representation);

}  // namespace internal_approximate_quantity

using internal_approximate_quantity::operator*;
using internal_approximate_quantity::operator""_⑴;
using internal_approximate_quantity::operator""_⑵;
using internal_approximate_quantity::operator""_⑶;
using internal_approximate_quantity::operator""_⑷;
using internal_approximate_quantity::operator""_⑸;
using internal_approximate_quantity::operator""_⑹;
using internal_approximate_quantity::operator""_⑺;
using internal_approximate_quantity::operator""_⑻;
using internal_approximate_quantity::operator""_⑼;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/approximate_quantity_body.hpp"
