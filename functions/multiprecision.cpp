#include "functions/multiprecision.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

// The designers of the Boost multiprecision library, in their confusion, have
// decided that nobody would use angles greater than about 1 / ε, so here I am,
// using a 2000-bit rational approximation of 2π to do a first round of
// reduction before calling their library.  Sigh.
cpp_rational const two_π(
    "85501617327051614281806889507174341018563881715359722726460240448778253362"
    "28033574412989010393845826881894399169689030062477829223652983133581682809"
    "07521562450510217951931573900693175198362946680952656259522439979797047648"
    "97394953111420464708282657298132250506331749517560155388580727329016311453"
    "71910",
    "13608005039951911862130979398099086168163224029757063948525789956741291980"
    "85326865353012588752613375780472488965315733291185793228273924312743154885"
    "30992776567452143477216307689797948122382208381130766696418513975898402329"
    "95420184160961245582170902313515516285224226723935852853290026779141303674"
    "08461");

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
