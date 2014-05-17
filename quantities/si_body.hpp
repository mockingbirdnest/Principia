#pragma once

namespace principia {
namespace si {
  template<typename D>
  inline quantities::Quantity<D> Yotta(quantities::Quantity<D> base) {
    return 1e24 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Zetta(quantities::Quantity<D> base) {
    return 1e21 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Exa(quantities::Quantity<D> base) {
    return 1e18 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Peta(quantities::Quantity<D> base) {
    return 1e15 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Tera(quantities::Quantity<D> base) {
    return 1e12 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Giga(quantities::Quantity<D> base) {
    return 1e9 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Mega(quantities::Quantity<D> base) {
    return 1e6 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Kilo(quantities::Quantity<D> base) {
    return 1e3 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Hecto(quantities::Quantity<D> base) {
    return 1e2 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Deca(quantities::Quantity<D> base) {
    return 1e1 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Deci(quantities::Quantity<D> base) {
    return 1e-1 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Centi(quantities::Quantity<D> base) {
    return 1e-2 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Milli(quantities::Quantity<D> base) {
    return 1e-3 * base;
  }
  template<typename D>
  quantities::Quantity<D> Micro(quantities::Quantity<D> base) {
    return 1e-6 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Nano(quantities::Quantity<D> base) {
    return 1e-9 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Pico(quantities::Quantity<D> base) {
    return 1e-12 * base;
  }
  template<typename D> 
  inline quantities::Quantity<D> Femto(quantities::Quantity<D> base) {
    return 1e-15 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Atto(quantities::Quantity<D> base) {
    return 1e-18 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Zepto(quantities::Quantity<D> base) {
    return 1e-21 * base;
  }
  template<typename D>
  inline quantities::Quantity<D> Yocto(quantities::Quantity<D> base) {
    return 1e-24 * base;
  }
}  // namespace si 
}  // namespace principia
