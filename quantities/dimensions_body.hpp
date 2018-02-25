
#pragma once

#include "quantities/dimensions.hpp"

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace internal_dimensions {

using base::not_constructible;

class ExponentSerializer : not_constructible {
 public:
  // Returns true if the exponent is in the range that we can serialize.
  static constexpr bool IsSerializable(std::int64_t exponent);

  // Returns the serialized representation of the exponent.  |position| is the
  // 0-based position of the dimension in the representation.
  static constexpr std::int64_t Representation(
      std::int64_t exponent,
      std::int64_t position);

 private:
  static constexpr std::int64_t min_exponent = -24;
  static constexpr std::int64_t max_exponent = 7;
  static constexpr std::int64_t exponent_mask = 0x1F;
  static constexpr std::int64_t exponent_bits = 5;
};

constexpr bool ExponentSerializer::IsSerializable(
    std::int64_t const exponent) {
  return exponent >= min_exponent && exponent <= max_exponent;
}

constexpr std::int64_t ExponentSerializer::Representation(
    std::int64_t const exponent,
    std::int64_t const position) {
  // For exponents in [-16, 7] this returns the representations
  // 0x10, 0x11, ... 0x00, ... 0x07.  For exponents in [-24, -17] this returns
  // the representations 0x08, 0x09, ... 0x0F.  The latter used to be reserved
  // for exponents in the range [8, 15] but we believe that we never used them,
  // and with polynomials in the monomial basis we need large negative
  // exponents.
  return (exponent & exponent_mask) << (position * exponent_bits);
}

template<std::int64_t LengthExponent,
         std::int64_t MassExponent,
         std::int64_t TimeExponent,
         std::int64_t CurrentExponent,
         std::int64_t TemperatureExponent,
         std::int64_t AmountExponent,
         std::int64_t LuminousIntensityExponent,
         std::int64_t AngleExponent>
struct Dimensions : not_constructible {
  enum {
    Length            = LengthExponent,
    Mass              = MassExponent,
    Time              = TimeExponent,
    Current           = CurrentExponent,
    Temperature       = TemperatureExponent,
    Amount            = AmountExponent,
    LuminousIntensity = LuminousIntensityExponent,
    Angle             = AngleExponent,
  };

  static std::int64_t constexpr representation =
      ExponentSerializer::Representation(LengthExponent, 0)            |
      ExponentSerializer::Representation(MassExponent, 1)              |
      ExponentSerializer::Representation(TimeExponent, 2)              |
      ExponentSerializer::Representation(CurrentExponent, 3)           |
      ExponentSerializer::Representation(TemperatureExponent, 4)       |
      ExponentSerializer::Representation(AmountExponent, 5)            |
      ExponentSerializer::Representation(LuminousIntensityExponent, 6) |
      ExponentSerializer::Representation(AngleExponent, 7);
};

template<typename Dimensions>
struct DimensionsAreSerializable {
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Length),
                "Invalid length exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Mass),
                "Invalid mass exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Time),
                "Invalid time exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Current),
                "Invalid current exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Temperature),
                "Invalid temperature exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Amount),
                "Invalid amount exponent");
  static_assert(ExponentSerializer::IsSerializable(
                    Dimensions::LuminousIntensity),
                "Invalid luminous intensity exponent");
  static_assert(ExponentSerializer::IsSerializable(Dimensions::Angle),
                "Invalid angle exponent");
};

template<typename Dimensions, int n>
struct DimensionsExponentiationGenerator : not_constructible {
  using Type =
      internal_dimensions::Dimensions<Dimensions::Length * n,
                                      Dimensions::Mass * n,
                                      Dimensions::Time * n,
                                      Dimensions::Current * n,
                                      Dimensions::Temperature * n,
                                      Dimensions::Amount * n,
                                      Dimensions::LuminousIntensity * n,
                                      Dimensions::Angle * n>;
};

template<typename Dimensions, int n>
struct DimensionsNthRootGenerator : not_constructible {
  static_assert((Dimensions::Length % n) == 0 &&
                (Dimensions::Mass % n) == 0 &&
                (Dimensions::Time % n) == 0 &&
                (Dimensions::Current % n) == 0 &&
                (Dimensions::Temperature % n) == 0 &&
                (Dimensions::Amount % n) == 0 &&
                (Dimensions::LuminousIntensity % n) == 0 &&
                (Dimensions::Angle % n) == 0,
                "Dimensions not suitable for Nth root");

  using Type =
      internal_dimensions::Dimensions<Dimensions::Length / n,
                                      Dimensions::Mass / n,
                                      Dimensions::Time / n,
                                      Dimensions::Current / n,
                                      Dimensions::Temperature / n,
                                      Dimensions::Amount / n,
                                      Dimensions::LuminousIntensity / n,
                                      Dimensions::Angle / n>;
};

template<typename LDimensions, typename RDimensions>
struct DimensionsProductGenerator : not_constructible {
  using Type = Dimensions<LDimensions::Length + RDimensions::Length,
                          LDimensions::Mass + RDimensions::Mass,
                          LDimensions::Time + RDimensions::Time,
                          LDimensions::Current + RDimensions::Current,
                          LDimensions::Temperature + RDimensions::Temperature,
                          LDimensions::Amount + RDimensions::Amount,
                          LDimensions::LuminousIntensity +
                              RDimensions::LuminousIntensity,
                          LDimensions::Angle + RDimensions::Angle>;
};

template<typename LDimensions, typename RDimensions>
struct DimensionsQuotientGenerator : not_constructible {
  using Type = Dimensions<LDimensions::Length - RDimensions::Length,
                          LDimensions::Mass - RDimensions::Mass,
                          LDimensions::Time - RDimensions::Time,
                          LDimensions::Current - RDimensions::Current,
                          LDimensions::Temperature - RDimensions::Temperature,
                          LDimensions::Amount - RDimensions::Amount,
                          LDimensions::LuminousIntensity -
                              RDimensions::LuminousIntensity,
                          LDimensions::Angle - RDimensions::Angle>;
};

}  // namespace internal_dimensions
}  // namespace quantities
}  // namespace principia
