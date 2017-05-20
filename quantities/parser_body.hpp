
#pragma once

#include "quantities/parser.hpp"

#include <array>
#include <string>

#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace internal_parser {

using si::AstronomicalUnit;
using si::Day;
using si::Degree;
using si::Kilo;
using si::Metre;
using si::Second;

template<typename Q>
struct ExtractDimensions {};

template<std::int64_t LengthExponent,
         std::int64_t MassExponent,
         std::int64_t TimeExponent,
         std::int64_t CurrentExponent,
         std::int64_t TemperatureExponent,
         std::int64_t AmountExponent,
         std::int64_t LuminousIntensityExponent,
         std::int64_t AngleExponent>
struct ExtractDimensions<
    Quantity<internal_quantities::Dimensions<LengthExponent,
                                             MassExponent,
                                             TimeExponent,
                                             CurrentExponent,
                                             TemperatureExponent,
                                             AmountExponent,
                                             LuminousIntensityExponent,
                                             AngleExponent>>> {
  constexpr static std::array<std::int64_t, 8> dimensions = {
      LengthExponent,
      MassExponent,
      TimeExponent,
      CurrentExponent,
      TemperatureExponent,
      AmountExponent,
      LuminousIntensityExponent,
      AngleExponent};
};

struct ParsedUnit {
  std::array<std::int64_t, 8> dimensions;
  double multiplier;
};

ParsedUnit ParseUnit(std::string const& s) {
  // Units of length.
  if (s == "m") {
    return {ExtractDimensions<Length>::dimensions, 1};
  } else if (s == "km") {
    return {ExtractDimensions<Length>::dimensions, Kilo(Metre) / Metre};
  } else if (s == "au") {
    return {ExtractDimensions<Length>::dimensions, AstronomicalUnit / Metre};
  // Units of time.
  } else if (s == "s") {
    return {ExtractDimensions<Time>::dimensions, 1};
  } else if (s == "d") {
    return {ExtractDimensions<Time>::dimensions, Day / Second};
  // Units of angle.
  } else if (s == "deg" || s == u8"°") {
    return {ExtractDimensions<Angle>::dimensions, Degree / Radian};
  } else if (s == "rad") {
    return {ExtractDimensions<Angle>::dimensions, 1};
  } else {
    LOG(FATAL) << "Unsupported unit " << s;
    base::noreturn();
  }
}

int ParseExponent(std::string const& s) {
  // Parse an int.
  char* interpreted_end;
  char const* const c_string = s.c_str();
  double const exponent = std::strtol(c_string, &interpreted_end, /*base=*/10);
  int const interpreted = interpreted_end - c_string;
  CHECK_LT(0, interpreted) << "invalid integer number " << s;
  return exponent;
}

ParsedUnit ParseExponentiationUnit(std::string const& s) {
  int const first_caret = s.find('^');
  if (first_caret == std::string::npos) {
    return ParseUnit(s);
  } else {
    int const first_nonblank = s.find_first_not_of(' ', last_slash + 1);
    CHECK_NE(std::string::npos, first_nonblank);
    int const last_nonblank = s.find_last_not_of(' ', last_slash - 1);
    CHECK_NE(std::string::npos, last_nonblank);
    auto const left = ParseUnit(s.substr(0, last_nonblank + 1));
    auto const right = ParseExponent(s.substr(first_nonblank));
  }
}

ParsedUnit ParseProductUnit(std::string const& s) {
  // For a product we are looking for a blank character that is not next to a
  // carret.
  int first_blank;
  int first_nonblank;
  int last_nonblank;
  for (;;) {
    first_blank = s.find(' ');
    if (first_blank == std::string::npos) {
      return ParseExponentiationUnit(s);
    } else {
      first_nonblank = s.find_first_not_of(' ', first_blank + 1);
      last_nonblank = s.find_last_not_of(' ', first_blank + 1);
      if ((first_nonblank == std::string::npos || s[first_nonblank] != '^') &&
          (last_nonblank == std::string::npos || s[last_nonblank] != '^')) {
        break;
      }
    }
  }
  auto const left = ParseProductUnit(s.substr(0, last_nonblank + 1));
  auto const right = ParseProductUnit(s.substr(first_nonblank));
}

ParsedUnit ParseQuotientUnit(std::string const& s) {
  // Look for the slash from the back to achieve proper associativity.
  int const last_slash = s.rfind('/');
  if (last_slash == std::string::npos) {
    // Not a quotient.
    return ParseProductUnit(s);
  } else {
    // A quotient.  Parse each half.
    int const first_nonblank = s.find_first_not_of(' ', last_slash + 1);
    CHECK_NE(std::string::npos, first_nonblank);
    int const last_nonblank = s.find_last_not_of(' ', last_slash - 1);
    CHECK_NE(std::string::npos, last_nonblank);
    auto const left = ParseProductUnit(s.substr(0, last_nonblank + 1));
    auto const right = ParseProductUnit(s.substr(first_nonblank));
  }
}

template<typename T>
T ParseQuantity(std::string const& s) {
  // Parse a double.
  char* interpreted_end;
  char const* const c_string = s.c_str();
  double const magnitude = std::strtod(c_string, &interpreted_end);
  int const interpreted = interpreted_end - c_string;
  CHECK_LT(0, interpreted) << "invalid floating-point number " << s;

  // Locate the unit.  It may be empty for a double.
  int const first_nonblank = s.find_first_not_of(' ', interpreted);
  int const last_nonblank = s.find_last_not_of(' ');
  std::string unit_string;
  if (first_nonblank != std::string::npos) {
    unit_string = s.substr(first_nonblank, last_nonblank - first_nonblank + 1);
  }

  ParsedUnit const unit = ParseQuotientUnit(unit_string);
  return magnitude * unit;
}

template<>
inline Length ParseUnit(std::string const& s) {
  if (s == "m") {
    return si::Metre;
  } else if (s == "km") {
    return si::Kilo(si::Metre);
  } else if (s == "au") {
    return si::AstronomicalUnit;
  } else {
    LOG(FATAL) << "Unsupported unit of length " << s;
    base::noreturn();
  }
}

template<>
inline Time ParseUnit(std::string const& s) {
  if (s == "s") {
    return si::Second;
  } else if (s == "d") {
    return si::Day;
  } else {
    LOG(FATAL) << "Unsupported unit of time " << s;
    base::noreturn();
  }
}

template<>
inline Angle ParseUnit(std::string const& s) {
  if (s == "deg" || s == u8"°") {
    return si::Degree;
  } else if (s == "rad") {
    return si::Radian;
  } else {
    LOG(FATAL) << "Unsupported unit of angle " << s;
    base::noreturn();
  }
}

template<>
inline AngularFrequency ParseUnit(std::string const& s) {
  return ParseQuotientUnit(s, &ParseUnit<Angle>, &ParseUnit<Time>);
}

template<>
inline Speed ParseUnit(std::string const& s) {
  return ParseQuotientUnit(s, &ParseUnit<Length>, &ParseUnit<Time>);
}

template<>
inline GravitationalParameter ParseUnit(std::string const& s) {
  return ParseQuotientUnit(s,
                           &ParseExponentiationUnit<Length, 3>,
                           &ParseExponentiationUnit<Time, 2>);
}

template<>
inline double ParseUnit<double>(std::string const& s) {
  CHECK(s.empty()) << s;
  return 1;
}

}  // namespace internal_parser
}  // namespace quantities
}  // namespace principia
