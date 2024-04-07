#include "functions/multiprecision.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

// The designers of the Boost multiprecision library, in their confusion, have
// decided that nobody would use angles greater than about 1 / ε, so here I am,
// using a 1000-bit rational approximation of 2π to do a first round of
// reduction before calling their library.  Sigh.
cpp_rational const two_π(
    "74778118356356393099550700487949553821548169522945304260092193625666650672"
    "12764554292012625975645707431952036396539900828099263122542425779444638301"
    "713",
    "11901307171524915725967286244795184945313912891969281885856459971598729125"
    "90781909108254956448232159392109804270301891472482302074906413637232243084"
    "548");

cpp_bin_float_50 Sin(cpp_rational const& α) {
  auto const k = static_cast<cpp_int>(α / two_π);
  auto const α_mod_2π = α - k * two_π;
  return sin(static_cast<cpp_bin_float_50>(α_mod_2π));
}

cpp_bin_float_50 Cos(cpp_rational const& α) {
  auto const k = static_cast<cpp_int>(α / two_π);
  auto const α_mod_2π = α - k * two_π;
  return cos(static_cast<cpp_bin_float_50>(α_mod_2π));
}

}  // namespace internal
}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
