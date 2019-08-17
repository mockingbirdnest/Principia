
#pragma once

#include "quantities/si.hpp"

#include <string>

namespace principia {
namespace quantities {
namespace si {

template<typename D>
std::string Format() {
  auto const format_unit = [](std::string const& name,
                              int const exponent) -> std::string {
    switch (exponent) {
      case 0:
        return "";
        break;
      case 1:
        return " " + name;
      default:
        return " " + name + "^" + std::to_string(exponent);
    }
  };

  // This string has a leading space if it's not empty.
  auto const format =
      format_unit("m", D::Length) + format_unit("kg", D::Mass) +
      format_unit("s", D::Time) + format_unit("A", D::Current) +
      format_unit("K", D::Temperature) + format_unit("mol", D::Amount) +
      format_unit("cd", D::LuminousIntensity) + format_unit("rad", D::Angle);

  if (format.empty()) {
    return format;
  } else {
    return format.substr(1, format.size() - 1);
  }
}

template<typename D>
constexpr Quantity<D> Yotta(Quantity<D> base) {
  return 1e24 * base;
}
template<typename D>
constexpr Quantity<D> Zetta(Quantity<D> base) {
  return 1e21 * base;
}
template<typename D>
constexpr Quantity<D> Exa(Quantity<D> base) {
  return 1e18 * base;
}
template<typename D>
constexpr Quantity<D> Peta(Quantity<D> base) {
  return 1e15 * base;
}
template<typename D>
constexpr Quantity<D> Tera(Quantity<D> base) {
  return 1e12 * base;
}
template<typename D>
constexpr Quantity<D> Giga(Quantity<D> base) {
  return 1e9 * base;
}
template<typename D>
constexpr Quantity<D> Mega(Quantity<D> base) {
  return 1e6 * base;
}
template<typename D>
constexpr Quantity<D> Kilo(Quantity<D> base) {
  return 1e3 * base;
}
template<typename D>
constexpr Quantity<D> Hecto(Quantity<D> base) {
  return 1e2 * base;
}
template<typename D>
constexpr Quantity<D> Deca(Quantity<D> base) {
  return 1e1 * base;
}
template<typename D>
constexpr Quantity<D> Deci(Quantity<D> base) {
  return 1e-1 * base;
}
template<typename D>
constexpr Quantity<D> Centi(Quantity<D> base) {
  return 1e-2 * base;
}
template<typename D>
constexpr Quantity<D> Milli(Quantity<D> base) {
  return 1e-3 * base;
}
template<typename D>
constexpr Quantity<D> Micro(Quantity<D> base) {
  return 1e-6 * base;
}
template<typename D>
constexpr Quantity<D> Nano(Quantity<D> base) {
  return 1e-9 * base;
}
template<typename D>
constexpr Quantity<D> Pico(Quantity<D> base) {
  return 1e-12 * base;
}
template<typename D>
constexpr Quantity<D> Femto(Quantity<D> base) {
  return 1e-15 * base;
}
template<typename D>
constexpr Quantity<D> Atto(Quantity<D> base) {
  return 1e-18 * base;
}
template<typename D>
constexpr Quantity<D> Zepto(Quantity<D> base) {
  return 1e-21 * base;
}
template<typename D>
constexpr Quantity<D> Yocto(Quantity<D> base) {
  return 1e-24 * base;
}

}  // namespace si
}  // namespace quantities
}  // namespace principia
