#pragma once

#include "quantities/parser.hpp"

#include <array>
#include <string>
#include <utility>

#include "quantities/astronomy.hpp"
#include "quantities/bipm.hpp"
#include "quantities/dimensions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace _parser {
namespace internal {

using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_bipm;
using namespace principia::quantities::_dimensions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

using RuntimeDimensions = std::array<std::int64_t, 8>;

template<typename Q>
struct ExtractDimensions {};

template<>
struct ExtractDimensions<double> {
  static constexpr RuntimeDimensions dimensions() {
    return {{0, 0, 0, 0, 0, 0, 0, 0}};
  }
};

template<std::int64_t LengthExponent,
         std::int64_t MassExponent,
         std::int64_t TimeExponent,
         std::int64_t CurrentExponent,
         std::int64_t TemperatureExponent,
         std::int64_t AmountExponent,
         std::int64_t LuminousIntensityExponent,
         std::int64_t AngleExponent>
struct ExtractDimensions<Quantity<Dimensions<LengthExponent,
                                             MassExponent,
                                             TimeExponent,
                                             CurrentExponent,
                                             TemperatureExponent,
                                             AmountExponent,
                                             LuminousIntensityExponent,
                                             AngleExponent>>> {
  static constexpr RuntimeDimensions dimensions() {
    return {{LengthExponent,
             MassExponent,
             TimeExponent,
             CurrentExponent,
             TemperatureExponent,
             AmountExponent,
             LuminousIntensityExponent,
             AngleExponent}};
  }
};

struct Unit {
  template<typename Q>
  explicit Unit(Q const& quantity);

  Unit(RuntimeDimensions&& dimensions, double scale);

  RuntimeDimensions dimensions;
  double scale;
};

template<typename Q>
Unit::Unit(Q const& quantity)
    : dimensions(ExtractDimensions<Q>::dimensions()),
      scale(quantity / si::Unit<Q>) {}

inline Unit::Unit(RuntimeDimensions&& dimensions, double const scale)
    : dimensions(dimensions),
      scale(scale) {}

inline Unit operator*(Unit const& left, Unit const& right) {
  RuntimeDimensions dimensions;
  for (std::int64_t i = 0; i < dimensions.size(); ++i) {
    dimensions[i] = left.dimensions[i] + right.dimensions[i];
  }
  return {std::move(dimensions), left.scale * right.scale};
}

inline Unit operator/(Unit const& left, Unit const& right) {
  RuntimeDimensions dimensions;
  for (std::int64_t i = 0; i < dimensions.size(); ++i) {
    dimensions[i] = left.dimensions[i] - right.dimensions[i];
  }
  return {std::move(dimensions), left.scale / right.scale};
}

inline Unit operator^(Unit const& left, int const exponent) {
  RuntimeDimensions dimensions;
  for (std::int64_t i = 0; i < dimensions.size(); ++i) {
    dimensions[i] = left.dimensions[i] * exponent;
  }
  return {std::move(dimensions), std::pow(left.scale, exponent)};
}

inline Unit ParseUnit(std::string const& s) {
  // Unitless quantities.
  if (s == "") {
    return Unit(1.0);
  // Units of length.
  } else if (s == "Å") {
    return Unit(Ångström);
  } else if (s == "μm") {
    return Unit(si::Micro(si::Metre));
  } else if (s == "mm") {
    return Unit(si::Milli(si::Metre));
  } else if (s == "cm") {
    return Unit(si::Centi(si::Metre));
  } else if (s == "m") {
    return Unit(si::Metre);
  } else if (s == "km") {
    return Unit(si::Kilo(si::Metre));
  } else if (s == "R🜨") {
    return Unit(TerrestrialEquatorialRadius);
  } else if (s == "R☉") {
    return Unit(SolarRadius);
  } else if (s == "au") {
    return Unit(AstronomicalUnit);
  // Units of mass.
  } else if (s == "kg") {
    return Unit(si::Kilogram);
  // Units of time.
  } else if (s == "ms") {
    return Unit(si::Milli(si::Second));
  } else if (s == "s") {
    return Unit(si::Second);
  } else if (s == "min") {
    return Unit(si::Minute);
  } else if (s == "h") {
    return Unit(si::Hour);
  } else if (s == "d") {
    return Unit(si::Day);
  // Units of gravitational parameter.
  } else if (s == "GM🜨") {
    return Unit(TerrestrialGravitationalParameter);
  } else if (s == "GM☉") {
    return Unit(SolarGravitationalParameter);
  // Units of power.
  } else if (s == "W") {
    return Unit(si::Watt);
  // Units of angle.
  } else if (s == "deg" || s == "°") {
    return Unit(si::Degree);
  } else if (s == "rad") {
    return Unit(si::Radian);
  // Units of solid angle.
  } else if (s == "sr") {
    return Unit(si::Steradian);
  } else {
    LOG(FATAL) << "Unsupported unit " << s;
    std::abort();
  }
}

inline int ParseExponent(std::string const& s) {
  // Parse an int.
  char* interpreted_end;
  char const* const c_string = s.c_str();
  int const exponent = std::strtol(c_string, &interpreted_end, /*base=*/10);
  std::int64_t const interpreted = interpreted_end - c_string;
  CHECK_LT(0, interpreted) << "invalid integer number " << s;
  return exponent;
}

inline Unit ParseExponentiationUnit(std::string const& s) {
  std::int64_t const first_caret = s.find('^');
  if (first_caret == std::string::npos) {
    return ParseUnit(s);
  } else {
    std::int64_t const first_nonblank =
        s.find_first_not_of(' ', first_caret + 1);
    CHECK_NE(std::string::npos, first_nonblank);
    std::int64_t const last_nonblank = s.find_last_not_of(' ', first_caret - 1);
    CHECK_NE(std::string::npos, last_nonblank);
    auto const left = ParseUnit(s.substr(0, last_nonblank + 1));
    auto const right = ParseExponent(s.substr(first_nonblank));
    return left ^ right;
  }
}

inline Unit ParseProductUnit(std::string const& s) {
  // For a product we are looking for a blank character that is not next to a
  // carret.
  std::int64_t first_blank;
  std::int64_t first_nonblank;
  std::int64_t last_nonblank;
  for (std::int64_t start = 0;; start = first_blank + 1) {
    first_blank = s.find(' ', start);
    if (first_blank == std::string::npos) {
      return ParseExponentiationUnit(s);
    } else {
      first_nonblank = s.find_first_not_of(' ', first_blank + 1);
      last_nonblank = s.find_last_not_of(' ', first_blank - 1);
      if ((first_nonblank == std::string::npos || s[first_nonblank] != '^') &&
          (last_nonblank == std::string::npos || s[last_nonblank] != '^')) {
        break;
      }
    }
  }
  auto const left = ParseExponentiationUnit(s.substr(0, last_nonblank + 1));
  auto const right = ParseProductUnit(s.substr(first_nonblank));
  return left * right;
}

inline Unit ParseQuotientUnit(std::string const& s) {
  // Look for the slash from the back to achieve proper associativity.
  std::int64_t const last_slash = s.rfind('/');
  if (last_slash == std::string::npos) {
    // Not a quotient.
    return ParseProductUnit(s);
  } else {
    // A quotient.  Parse each half.  Note that there may not be a left half for
    // input like 1.23 / s.
    std::int64_t const first_nonblank =
        s.find_first_not_of(' ', last_slash + 1);
    CHECK_NE(std::string::npos, first_nonblank);
    std::size_t const last_nonblank =
        last_slash == 0 ? std::string::npos
                        : s.find_last_not_of(' ', last_slash - 1);
    auto const left = last_nonblank == std::string::npos
                          ? Unit(1.0)
                          : ParseQuotientUnit(s.substr(0, last_nonblank + 1));
    auto const right = ParseExponentiationUnit(s.substr(first_nonblank));
    return left / right;
  }
}

template<typename Q>
Q ParseQuantity(std::string const& s) {
  // Parse a double.
  char* interpreted_end;
  char const* const c_string = s.c_str();
  double const magnitude = std::strtod(c_string, &interpreted_end);
  std::int64_t const interpreted = interpreted_end - c_string;
  CHECK_LT(0, interpreted) << "invalid floating-point number " << s;

  // Locate the unit.  It may be empty for a double.
  std::int64_t const first_nonblank = s.find_first_not_of(' ', interpreted);
  std::int64_t const last_nonblank = s.find_last_not_of(' ');
  std::string unit_string;
  if (first_nonblank != std::string::npos) {
    unit_string = s.substr(first_nonblank, last_nonblank - first_nonblank + 1);
  }

  Unit const unit = ParseQuotientUnit(unit_string);
  CHECK(ExtractDimensions<Q>::dimensions() == unit.dimensions) << unit_string;
  return magnitude * unit.scale * si::Unit<Q>;
}

}  // namespace internal
}  // namespace _parser
}  // namespace quantities
}  // namespace principia
