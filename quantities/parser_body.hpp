#pragma once

#include "quantities/parser.hpp"

#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

namespace {

// The patterns that we parse here have the form:
//    U^n/V^m
//  when U and V are (final) unit names and n and m are integers.

template<typename T>
using ParseUnitFunction = T(*)(std::string const& s);

template<typename T, int exponent>
Exponentiation<T, exponent> ParseExponentiationUnit(std::string const& s) {
  int const first_carret = s.find('^');
  int const last_nonblank = s.find_last_not_of(' ', first_carret - 1);
  CHECK_NE(std::string::npos, last_nonblank);

  char* interpreted_end;
  const char* interpreted_begin = s.c_str() + first_carret + 1;
  int const actual_exponent = std::strtol(interpreted_begin,
                                          &interpreted_end,
                                          10);
  int const interpreted = interpreted_end - interpreted_begin;
  CHECK_LT(0, interpreted) << "invalid exponent " << s;
  CHECK_EQ(exponent, actual_exponent);

  return Pow<exponent>(ParseUnit<T>(s.substr(0, last_nonblank + 1)));
}

template<typename TNumerator, typename TDenominator>
Quotient<TNumerator, TDenominator> ParseQuotientUnit(
    std::string const& s,
    ParseUnitFunction<TNumerator> parse_numerator_unit,
    ParseUnitFunction<TDenominator> parse_denominator_unit) {
  int const first_slash = s.find('/');
  int const first_nonblank = s.find_first_not_of(' ', first_slash + 1);
  CHECK_NE(std::string::npos, first_nonblank);
  int const last_nonblank = s.find_last_not_of(' ', first_slash - 1);
  CHECK_NE(std::string::npos, last_nonblank);
  return parse_numerator_unit(s.substr(0, last_nonblank + 1)) /
         parse_denominator_unit(s.substr(first_nonblank));
}

}  // namespace

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

  T const unit = ParseUnit<T>(unit_string);
  return magnitude * unit;
}

template<>
Length ParseUnit(std::string const& s) {
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
Time ParseUnit(std::string const& s) {
  if (s == "s") {
    return si::Second;
  } else if (s == "d") {
    return si::Day;
  } else {
    LOG(FATAL) << "Unsupported unit of speed " << s;
    base::noreturn();
  }
}

template<>
Angle ParseUnit(std::string const& s) {
  if (s == "deg" || s == "°") {
    return si::Degree;
  } else if (s == "rad") {
    return si::Radian;
  } else {
    LOG(FATAL) << "Unsupported unit of angle " << s;
    base::noreturn();
  }
}

template<>
Speed ParseUnit(std::string const& s) {
  return ParseQuotientUnit(s, &ParseUnit<Length>, &ParseUnit<Time>);
}

template<>
GravitationalParameter ParseUnit(std::string const& s) {
  return ParseQuotientUnit(s,
                           &ParseExponentiationUnit<Length, 3>,
                           &ParseExponentiationUnit<Time, 2>);
}

template<>
double ParseUnit<double>(std::string const& s) {
  CHECK(s.empty()) << s;
  return 1;
}

}  // namespace quantities
}  // namespace principia
