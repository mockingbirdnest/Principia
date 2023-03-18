#include "integrators/adams_moulton_integrator.hpp"

namespace principia {
namespace integrators {
namespace _adams_moulton_integrator {
namespace internal {

// For the Adams-Moulton integrators, see
// http://www.scholarpedia.org/article/Linear_multistep_method#Methods_for_Hamiltonian_systems  // NOLINT
// and
// https://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods.  // NOLINT
// The latter has a general formula for the coefficients and was used to
// compute the high-order integrators.

template<>
inline AdamsMoulton<1> const& AdamsMoultonOrder<1>() {
  static AdamsMoulton<1> const integrator{{1.0}, 1.0};
  return integrator;
}

template<>
inline AdamsMoulton<2> const& AdamsMoultonOrder<2>() {
  static AdamsMoulton<2> const integrator{{1.0, 1.0}, 2.0};
  return integrator;
}

template<>
inline AdamsMoulton<3> const& AdamsMoultonOrder<3>() {
  static AdamsMoulton<3> const integrator{{5.0, 8.0, -1.0}, 12.0};
  return integrator;
}

template<>
inline AdamsMoulton<4> const& AdamsMoultonOrder<4>() {
  static AdamsMoulton<4> const integrator{{9.0, 19.0, -5.0, 1.0}, 24.0};
  return integrator;
}

template<>
inline AdamsMoulton<5> const& AdamsMoultonOrder<5>() {
  static AdamsMoulton<5> const integrator{{251.0, 646.0, -264.0, 106.0, -19.0},
                                          720.0};
  return integrator;
}

template<>
inline AdamsMoulton<6> const& AdamsMoultonOrder<6>() {
  static AdamsMoulton<6> const integrator{
      {475.0, 1427.0, -798.0, 482.0, -173.0, 27.0}, 1440.0};
  return integrator;
}

template<>
inline AdamsMoulton<7> const& AdamsMoultonOrder<7>() {
  static AdamsMoulton<7> const integrator{
      {19087.0, 65112.0, -46461.0, 37504.0, -20211.0, 6312.0, -863.0}, 60480.0};
  return integrator;
}

template<>
inline AdamsMoulton<8> const& AdamsMoultonOrder<8>() {
  static AdamsMoulton<8> const integrator{{36799.0,
                                           139849.0,
                                           -121797.0,
                                           123133.0,
                                           -88547.0,
                                           41499.0,
                                           -11351.0,
                                           1375.0},
                                          120960.0};
  return integrator;
}

template<>
inline AdamsMoulton<9> const& AdamsMoultonOrder<9>() {
  static AdamsMoulton<9> const integrator{{1070017.0,
                                           4467094.0,
                                           -4604594.0,
                                           5595358.0,
                                           -5033120.0,
                                           3146338.0,
                                           -1291214.0,
                                           312874.0,
                                           -33953.0},
                                          3628800.0};
  return integrator;
}

template<>
inline AdamsMoulton<10> const& AdamsMoultonOrder<10>() {
  static AdamsMoulton<10> const integrator{{2082753.0,
                                            9449717.0,
                                            -11271304.0,
                                            16002320.0,
                                            -17283646.0,
                                            13510082.0,
                                            -7394032.0,
                                            2687864.0,
                                            -583435.0,
                                            57281.0},
                                           7257600.0};
  return integrator;
}

template<>
inline AdamsMoulton<11> const& AdamsMoultonOrder<11>() {
  static AdamsMoulton<11> const integrator{{134211265.0,
                                            656185652.0,
                                            -890175549.0,
                                            1446205080.0,
                                            -1823311566.0,
                                            1710774528.0,
                                            -1170597042.0,
                                            567450984.0,
                                            -184776195.0,
                                            36284876.0,
                                            -3250433.0},
                                           479001600.0};
  return integrator;
}

template<>
inline AdamsMoulton<12> const& AdamsMoultonOrder<12>() {
  static AdamsMoulton<12> const integrator{{262747265.0,
                                            1374799219.0,
                                            -2092490673.0,
                                            3828828885.0,
                                            -5519460582.0,
                                            6043521486.0,
                                            -4963166514.0,
                                            3007739418.0,
                                            -1305971115.0,
                                            384709327.0,
                                            -68928781.0,
                                            5675265.0},
                                           958003200.0};
  return integrator;
}

template<>
inline AdamsMoulton<13> const& AdamsMoultonOrder<13>() {
  static AdamsMoulton<13> const integrator{{703604254357.0,
                                            3917551216986.0,
                                            -6616420957428.0,
                                            13465774256510.0,
                                            -21847538039895.0,
                                            27345870698436.0,
                                            -26204344465152.0,
                                            19058185652796.0,
                                            -10344711794985.0,
                                            4063327863170.0,
                                            -1092096992268.0,
                                            179842822566.0,
                                            -13695779093.0},
                                           2615348736000.0};
  return integrator;
}

template<>
inline AdamsMoulton<14> const& AdamsMoultonOrder<14>() {
  static AdamsMoulton<14> const integrator{{1382741929621.0,
                                            8153167962181.0,
                                            -15141235084110.0,
                                            33928990133618.0,
                                            -61188680131285.0,
                                            86180228689563.0,
                                            -94393338653892.0,
                                            80101021029180.0,
                                            -52177910882661.0,
                                            25620259777835.0,
                                            -9181635605134.0,
                                            2268078814386.0,
                                            -345457086395.0,
                                            24466579093.0},
                                           5230697472000.0};
  return integrator;
}

}  // namespace internal
}  // namespace _adams_moulton_integrator
}  // namespace integrators
}  // namespace principia
