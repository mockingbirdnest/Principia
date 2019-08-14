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
class ApproximateQuantity {
 public:
  Quantity min() const;
  Quantity max() const;

 private:
  ApproximateQuantity(std::string const& representation,
                      Quantity const& unit,
                      double min_multiplier,
                      double max_multiplier);

  // The original representation.
  std::string representation_;

  // The unit for the approximate quantity.  Could for instance be Degree for an
  // angle.
  Quantity unit_;

  // The interval for the approximate quantity, expressed as multiples of unit_.
  double min_multiplier_;
  double max_multiplier_;

  template<typename Left, typename Right>
  friend Product<Left, Right> operator*(ApproximateQuantity<Left> const& left,
                                        Right const& right);

  friend ApproximateQuantity<double> operator""_⑴(char const* representation);
};

template<typename Left, typename Right>
Product<Left, Right> operator*(ApproximateQuantity<Left> const& left,
                               Right const& right);

ApproximateQuantity<double> operator""_⑴(char const* representation);

}  // namespace internal_approximate_quantity
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/approximate_quantity_body.hpp"
