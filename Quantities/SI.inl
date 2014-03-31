// SIUnits.ipp

#pragma once

namespace Principia {
namespace SI {
  template<typename D>
  inline Quantities::Quantity<D> Yotta(Quantities::Quantity<D> base) {
    return 1e24 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Zetta(Quantities::Quantity<D> base) {
    return 1e21 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Exa(Quantities::Quantity<D> base) {
    return 1e18 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Peta(Quantities::Quantity<D> base) {
    return 1e15 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Tera(Quantities::Quantity<D> base) {
    return 1e12 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Giga(Quantities::Quantity<D> base) {
    return 1e9 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Mega(Quantities::Quantity<D> base) {
    return 1e6 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Kilo(Quantities::Quantity<D> base) {
    return 1e3 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Hecto(Quantities::Quantity<D> base) {
    return 1e2 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Deca(Quantities::Quantity<D> base) {
    return 1e1 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Deci(Quantities::Quantity<D> base) {
    return 1e-1 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Centi(Quantities::Quantity<D> base) {
    return 1e-2 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Milli(Quantities::Quantity<D> base) {
    return 1e-3 * base;
  }
  template<typename D> Quantities::Quantity<D> Micro(Quantities::Quantity<D> base) {
    return 1e-6 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Nano(Quantities::Quantity<D> base) {
    return 1e-9 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Pico(Quantities::Quantity<D> base) {
    return 1e-12 * base;
  }
  template<typename D> 
  inline Quantities::Quantity<D> Femto(Quantities::Quantity<D> base) {
    return 1e-15 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Atto(Quantities::Quantity<D> base) {
    return 1e-18 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Zepto(Quantities::Quantity<D> base) {
    return 1e-21 * base;
  }
  template<typename D>
  inline Quantities::Quantity<D> Yocto(Quantities::Quantity<D> base) {
    return 1e-24 * base;
  }
}
}
