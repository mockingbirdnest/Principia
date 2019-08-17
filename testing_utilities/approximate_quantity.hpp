#pragma once

#include <string>
#include <string_view>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace internal_approximate_quantity {

using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;

template<typename Quantity>
class ApproximateQuantity;

template<typename Dimensions>
class ApproximateQuantity<Quantity<Dimensions>> {
 public:
  Quantity<Dimensions> min() const;
  Quantity<Dimensions> max() const;

  std::string DebugString() const;

 private:
  ApproximateQuantity(std::string const& representation,
                      int ulp,
                      bool negated,
                      double min_multiplier,
                      double max_multiplier,
                      Quantity<Dimensions> const& unit);

  // The original representation.
  std::string representation_;
  int ulp_;
  bool negated_;

  // The interval for the approximate quantity, expressed as multiples of unit_.
  double min_multiplier_;
  double max_multiplier_;

  // The unit for the approximate quantity.  Could for instance be Degree for an
  // angle.
  Quantity<Dimensions> unit_;

  template<typename Right>
  friend ApproximateQuantity<Right> operator-(
      ApproximateQuantity<Right> const& right);
  template<typename Left, typename RDimensions>
  friend ApproximateQuantity<Product<Left, Quantity<RDimensions>>> operator*(
      ApproximateQuantity<Left> const& left,
      Quantity<RDimensions> const& right);
  template<typename Left, typename RDimensions>
  friend ApproximateQuantity<Quotient<Left, Quantity<RDimensions>>> operator/(
      ApproximateQuantity<Left> const& left,
      Quantity<RDimensions> const& right);
};

template<>
class ApproximateQuantity<double> {
 public:
  static ApproximateQuantity<double> Parse(std::string_view representation,
                                           int ulp);

  double min() const;
  double max() const;

  std::string DebugString() const;

 private:
  ApproximateQuantity(std::string_view representation,
                      int ulp,
                      bool negated,
                      double min_multiplier,
                      double max_multiplier);

  // The original representation.
  std::string representation_;
  int ulp_;
  bool negated_;

  // The interval for the approximate quantity.
  double min_multiplier_;
  double max_multiplier_;

  static constexpr double unit_ = 1;

  template<typename Right>
  friend ApproximateQuantity<Right> operator-(
      ApproximateQuantity<Right> const& right);
  template<typename Left, typename RDimensions>
  friend ApproximateQuantity<Product<Left, Quantity<RDimensions>>> operator*(
      ApproximateQuantity<Left> const& left,
      Quantity<RDimensions> const& right);
  template<typename Left, typename RDimensions>
  friend ApproximateQuantity<Quotient<Left, Quantity<RDimensions>>> operator/(
      ApproximateQuantity<Left> const& left,
      Quantity<RDimensions> const& right);
};

template<typename Right>
ApproximateQuantity<Right> operator+(ApproximateQuantity<Right> const& right);
template<typename Right>
ApproximateQuantity<Right> operator-(ApproximateQuantity<Right> const& right);
template<>
ApproximateQuantity<double> operator-(ApproximateQuantity<double> const& right);
template<typename Left, typename RDimensions>
ApproximateQuantity<Product<Left, Quantity<RDimensions>>> operator*(
    ApproximateQuantity<Left> const& left,
    Quantity<RDimensions> const& right);
template<typename Left, typename RDimensions>
ApproximateQuantity<Quotient<Left, Quantity<RDimensions>>> operator/(
    ApproximateQuantity<Left> const& left,
    Quantity<RDimensions> const& right);

template<typename Quantity>
std::ostream& operator<<(std::ostream& out,
                         ApproximateQuantity<Quantity> const& q);

// The 🄐 to 🄕 operators are only for hexadecimal literals.
ApproximateQuantity<double> operator""_⑴(char const* representation);
ApproximateQuantity<double> operator""_⑵(char const* representation);
ApproximateQuantity<double> operator""_⑶(char const* representation);
ApproximateQuantity<double> operator""_⑷(char const* representation);
ApproximateQuantity<double> operator""_⑸(char const* representation);
ApproximateQuantity<double> operator""_⑹(char const* representation);
ApproximateQuantity<double> operator""_⑺(char const* representation);
ApproximateQuantity<double> operator""_⑻(char const* representation);
ApproximateQuantity<double> operator""_⑼(char const* representation);
ApproximateQuantity<double> operator""_🄐(char const* representation);
ApproximateQuantity<double> operator""_🄑(char const* representation);
ApproximateQuantity<double> operator""_🄒(char const* representation);
ApproximateQuantity<double> operator""_🄓(char const* representation);
ApproximateQuantity<double> operator""_🄔(char const* representation);
ApproximateQuantity<double> operator""_🄕(char const* representation);

}  // namespace internal_approximate_quantity

using internal_approximate_quantity::ApproximateQuantity;
using internal_approximate_quantity::operator*;
using internal_approximate_quantity::operator/;
using internal_approximate_quantity::operator""_⑴;
using internal_approximate_quantity::operator""_⑵;
using internal_approximate_quantity::operator""_⑶;
using internal_approximate_quantity::operator""_⑷;
using internal_approximate_quantity::operator""_⑸;
using internal_approximate_quantity::operator""_⑹;
using internal_approximate_quantity::operator""_⑺;
using internal_approximate_quantity::operator""_⑻;
using internal_approximate_quantity::operator""_⑼;
using internal_approximate_quantity::operator""_🄐;
using internal_approximate_quantity::operator""_🄑;
using internal_approximate_quantity::operator""_🄒;
using internal_approximate_quantity::operator""_🄓;
using internal_approximate_quantity::operator""_🄔;
using internal_approximate_quantity::operator""_🄕;

}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/approximate_quantity_body.hpp"
