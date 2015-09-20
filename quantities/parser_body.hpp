#pragma once

#include "quantities/parser.hpp"

#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

template<typename Q>
Q ParseQuantity(std::string const& s) {
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

  Q const unit = ParseUnit<Q>(unit_string);
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
Speed ParseUnit(std::string const& s) {
  if (s == "m/s") {
    return si::Metre / si::Second;
  } else if (s == "km/s") {
    return si::Kilo(si::Metre) / si::Second;
  } else if (s == "km/d") {
    return si::Kilo(si::Metre) / si::Day;
  } else if (s == "au/d") {
    return si::AstronomicalUnit / si::Day;
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
GravitationalParameter ParseUnit(std::string const& s) {
  if (s == "m^3/s^2") {
    return Pow<3>(si::Metre) / Pow<2>(si::Second);
  } else if (s == "km^3/s^2") {
    return Pow<3>(si::Kilo(si::Metre)) / Pow<2>(si::Second);
  } else if (s == "km^3/d^2") {
    return Pow<3>(si::Kilo(si::Metre)) / Pow<2>(si::Day);
  } else if (s == "au^3/d^2") {
    return Pow<3>(si::AstronomicalUnit) / Pow<2>(si::Day);
  } else {
    LOG(FATAL) << "Unsupported unit of gravitational parameter " << s;
    base::noreturn();
  }
}

template<>
double ParseUnit<double>(std::string const& s) {
  CHECK(s.empty()) << s;
  return 1;
}

}  // namespace quantities
}  // namespace principia
