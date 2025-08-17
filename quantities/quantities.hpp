#pragma once

#include <pmmintrin.h>

#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

#include "base/macros.hpp"  // 🧙 For CONSTEXPR_NAN.
#include "base/not_constructible.hpp"
#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "boost/multiprecision/number.hpp"
#include "quantities/dimensions.hpp"
#include "quantities/generators.hpp"
#include "serialization/quantities.pb.h"

namespace principia {
namespace quantities {
namespace _quantities {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_tags;
using namespace principia::quantities::_dimensions;
using namespace principia::quantities::_generators;

template<typename D>
class Quantity;

// Base quantities
using Length            = Quantity<Dimensions<1, 0, 0, 0, 0, 0, 0, 0>>;
using Mass              = Quantity<Dimensions<0, 1, 0, 0, 0, 0, 0, 0>>;
using Time              = Quantity<Dimensions<0, 0, 1, 0, 0, 0, 0, 0>>;
using Current           = Quantity<Dimensions<0, 0, 0, 1, 0, 0, 0, 0>>;
using Temperature       = Quantity<Dimensions<0, 0, 0, 0, 1, 0, 0, 0>>;
using Amount            = Quantity<Dimensions<0, 0, 0, 0, 0, 1, 0, 0>>;
using LuminousIntensity = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 1, 0>>;
// We strongly type angles.
using Angle             = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 1>>;

// `Product` and `Quotient` are not exported from this namespace.  Instead they
// are defined as the result types of `operator*` and `operator/`.
template<typename Left, typename Right>
using Product = typename ProductGenerator<Left, Right>::Type;
template<typename Left, typename Right>
using Quotient = typename QuotientGenerator<Left, Right>::Type;

template<typename D>
class Quantity final {
 public:
  using Dimensions = D;

  constexpr Quantity() = default;
  explicit constexpr Quantity(uninitialized_t);

  constexpr friend auto operator<=>(Quantity const& left,
                                    Quantity const& right) = default;

  constexpr Quantity operator+() const;
  constexpr Quantity operator-() const;
  constexpr Quantity operator+(Quantity const& right) const;
  constexpr Quantity operator-(Quantity const& right) const;

  constexpr Quantity operator*(double right) const;
  constexpr Quantity operator/(double right) const;

  Quantity& operator+=(Quantity const& right);
  Quantity& operator-=(Quantity const& right);
  Quantity& operator*=(double right);
  Quantity& operator/=(double right);

  void WriteToMessage(not_null<serialization::Quantity*> message) const;
  static Quantity ReadFromMessage(serialization::Quantity const& message);

 private:
  explicit constexpr Quantity(double magnitude);
  double magnitude_ = 0;

  template<typename LDimensions, typename RDimensions>
  friend constexpr Product<Quantity<LDimensions>,
                           Quantity<RDimensions>> operator*(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename LDimensions, typename RDimensions>
  friend constexpr Quotient<Quantity<LDimensions>,
                            Quantity<RDimensions>> operator/(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr Quantity<RDimensions> operator*(
      double left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr Quotient<double, Quantity<RDimensions>> operator/(
      double left,
      Quantity<RDimensions> const& right);

  template<typename Q>
  friend constexpr Q SIUnit();

  template<typename Dimensions>
  friend __m128d ToM128D(Quantity<Dimensions> x);
};

template<typename LDimensions, typename RDimensions>
constexpr Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename LDimensions, typename RDimensions>
constexpr Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
operator/(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quantity<RDimensions>
operator*(double, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quotient<double, Quantity<RDimensions>>
operator/(double, Quantity<RDimensions> const&);

// Used for implementing `si::Unit`.  Don't call directly, don't export from
// this namespace.  Defined here to break circular dependencies.
template<typename Q>
constexpr Q SIUnit() { return Q(1); };

inline __m128d ToM128D(double x);

template<typename Dimensions>
__m128d ToM128D(Quantity<Dimensions> x);

// A positive infinity of `Q`.
template<typename Q>
constexpr Q Infinity = SIUnit<Q>() * std::numeric_limits<double>::infinity();
// A quiet NaN of `Q`.
template <typename Q>
CONSTEXPR_NAN Q NaN = SIUnit<Q>() * std::numeric_limits<double>::quiet_NaN();

template<typename Q>
constexpr bool IsFinite(Q const& x);


template<typename D>
std::string Format();

std::string DebugString(
    double number,
    int precision = std::numeric_limits<double>::max_digits10);
template<typename N>
  requires is_number<N>::value
std::string DebugString(
    N const& number,
    int precision = std::numeric_limits<double>::max_digits10);
template<typename D>
std::string DebugString(
    Quantity<D> const& quantity,
    int precision = std::numeric_limits<double>::max_digits10);

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity);

}  // namespace internal

using internal::Amount;
using internal::Angle;
using internal::Current;
using internal::DebugString;
using internal::Format;
using internal::Infinity;
using internal::IsFinite;
using internal::Length;
using internal::LuminousIntensity;
using internal::Mass;
using internal::NaN;
using internal::Quantity;
using internal::Temperature;
using internal::Time;
using internal::ToM128D;

}  // namespace _quantities
}  // namespace quantities
}  // namespace principia

#include "quantities/quantities_body.hpp"
