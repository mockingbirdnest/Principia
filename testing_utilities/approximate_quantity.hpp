#pragma once

#include <string>
#include <string_view>

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace testing_utilities {
namespace _approximate_quantity {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Quantity>
class ApproximateQuantity;

template<typename Dimensions>
class ApproximateQuantity<Quantity<Dimensions>> {
 public:
  Quantity<Dimensions> min() const;
  Quantity<Dimensions> max() const;

  Quantity<Dimensions> unit() const;
  bool has_trivial_unit() const;

  // Returns the distance of the argument from the centre of the interval for
  // this object, in ulps.
  double UlpDistance(Quantity<Dimensions> const& q) const;

  std::string DebugString() const;

 private:
  ApproximateQuantity(std::string representation,
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

  // Returns the distance of the argument from the centre of the interval for
  // this object, in ulps.
  double UlpDistance(double d) const;

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

  double const unit_ = 1;

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

// The internal struct |NumericLiteral| and the operator""_ are to be used in
// concert to allow the syntax 1.234_(1), or, in principle, 1.234_(123) for
// approximate quantities.
struct NumericLiteral {
  ApproximateQuantity<double> operator()(int);
  char const* representation;
};

NumericLiteral operator""_(char const* representation);

}  // namespace internal

using internal::ApproximateQuantity;
using internal::operator*;
using internal::operator/;
using internal::operator""_;

}  // namespace _approximate_quantity
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/approximate_quantity_body.hpp"
