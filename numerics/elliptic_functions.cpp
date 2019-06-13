
#include "glog/logging.h"

namespace principia {
namespace numerics {

// Double precision subroutine to compute three Jacobian elliptic functions
// simultaneously
//
//   For general argument: -infty < u < infty
//
//     Reference: T. Fukushima, (2012) Numer. Math.
//     DOI 10.1007/s00211-012-0498-0
//       "Precise and Fast Computation of Jacobian Elliptic Functions by
//        Conditional Duplication"
//
//     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
//
//     Used subprograms: scd2, elk
//
//     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
//
//     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
//
void Gscd(double const u, double const mc, double& s, double& c, double& d) {
  double m, kc, ux, k, kh, kh3, kh5, kh7, k2, k3, k4, sx, cx, dx;
  double elk;

  m = 1.0 - mc;
  kc = sqrt(mc);
  ux = abs(u);
  if (ux < 0.785) {
    Scd2(ux, mc, s, c, d)
  } else {
    k = elk(mc);
    kh = k * 0.5;
    kh3 = k * 1.5;
    kh5 = k * 2.5;
    kh7 = k * 3.5;
    k2 = k * 2.0;
    k3 = k * 3.0;
    k4 = k * 4.0;
    ux = ux - k4 * static_cast<double>(static_cast<int>(ux / k4));
    if (ux < kh) {
      Scd2(ux, mc, s, c, d);
    } else if (ux < k) {
      ux = k - ux;
      Scd2(ux, mc, s, c, d);
      sx = c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh3) {
      ux = ux - k;
      Scd2(ux, mc, s, c, d);
      sx = c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < k2) {
      ux = k2 - ux;
      Scd2(ux, mc, s, c, d);
      c = -c;
    } else if (ux < kh5) {
      ux = ux - k2;
      Scd2(ux, mc, s, c, d);
      s = -s;
      c = -c;
    } else if (ux < k3) {
      ux = k3 - ux;
      Scd2(ux, mc, s, c, d);
      sx = -c / d;
      c = -kc * s / d;
      s = sx;
      d = kc / d;
    } else if (ux < kh7) {
      ux = ux - k3;
      Scd2(ux, mc, s, c, d);
      sx = -c / d;
      c = kc * s / d;
      s = sx;
      d = kc / d;
    } else {
      ux = k4 - ux;
      Scd2(ux, mc, s, c, d);
      s = -s;
    }
  }
  if (u < 0.0) {
    s = -s
  }
}

// Double precision subroutine to compute three Jacobian elliptic functions simultaneously
//
//   For limited argument: 0 <= u < K/2
//
//     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
//       "Precise and Fast Computation of Jacobian Elliptic Functions by
//        Conditional Duplication"
//
//     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
//
//     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
//
//     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
//
void Scd2(double const u, double const mc, double& s, double& c, double& d) {
  double B10, B11, B20, B21, B22, m, uA, uT, u0, v, a, b, y, z, my, mc2, m2, xz,
      w;
  int n, j, i;
  constexpr double B10 = 1.0 / 24.0;
  constexpr double B11 = 1.0 / 6.0;
  constexpr double B20 = 1.0 / 720.0;
  constexpr double B21 = 11.0 / 180.0;
  constexpr double B22 = 1.0 / 45.0;
  m = 1.0 - mc;
  uA = 1.76269 + mc * 1.16357;
  uT = 5.217e-3 - m * 2.143e-3;
  u0 = u;
  for (int n = 0; n <= 20; ++n) {
    if (u0 < uT) {
      break;
    }
    LOG_IF(FATAL, n == 20) << "(scd2) Too large input argument: u=" << u;
    u0 = u0 * 0.5;
  }
  v = u0 * u0;
  a = 1.0;
  b = v * (0.5 - v * (B10 + m * B11 - v * (B20 + m * (B21 + m * B22))));
  if (u < uA) {
    for (int j = 1; j <= n; ++j) {
      y = b * (a * 2.0 - b);
      z = a * a;
      my = m * y;
      b = (y * 2.0) * (z - my);
      a = z * z - my * y;
    }
  } else {
    for (int j = 1; j <= n; ++j) {
      y = b * (a * 2.0 - b);
      z = a * a;
      my = m * y;
      if (z < my * 2.0) {
        c = a - b;
        mc2 = mc * 2.0;
        m2 = m * 2.0;
        for (int i = j; i <= n; ++i) {
          x = c * c;
          z = a * a;
          w = m * x * x - mc * z * z;
          xz = x * z;
          c = mc2 * xz + w;
          a = m2 * xz - w;
        }
        c = c / a;
        x = c * c;
        s = sqrt(1.0 - x);
        d = sqrt(mc + m * x);
        return;
      }
      b = (y * 2.0) * (z - my);
      a = z * z - my * y
    }
  }
  b = b / a;
  y = b * (2.0 - b);
  c = 1.0 - b;
  s = sqrt(y);
  d = sqrt(1.0 - m * y);
}

//  Double precision complete elliptic integral of the first kind
//
//     Reference: T. Fukushima, (2009) Celest. Mech. Dyn. Astron. 105, 305-328
//        "Fast Computation of Complete Elliptic Integrlals and Jacobian
//         Elliptic Functions"
//
//     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
//
//     Inputs: mc   = complementary parameter 0 <= mc   <= 1
//
//     Output: elk
//
double Elk(double const mc) {
  double m, mx, P, Q;
  double kkc, nome;

  constexpr double D1 = 1.0 / 16.0;
  constexpr double D2 = 1.0 / 32.0;
  constexpr double D3 = 21.0 / 1024.0;
  constexpr double D4 = 31.0 / 2048.0;
  constexpr double D5 = 6257.0 / 524288.0;
  constexpr double D6 = 10293.0 / 1048576.0;
  constexpr double D7 = 279025.0 / 33554432.0;
  constexpr double D8 = 483127.0 / 67108864.0;
  constexpr double D9 = 435506703.0 / 68719476736.0;
  constexpr double D10 = 776957575.0 / 137438953472.0;
  constexpr double D11 = 22417045555.0 / 4398046511104.0;
  constexpr double D12 = 40784671953.0 / 8796093022208.0;
  constexpr double D13 = 9569130097211.0 / 2251799813685248.0;
  constexpr double D14 = 17652604545791.0 / 4503599627370496.0;

  static bool first = true;
  static double mcold, PIHALF, PIINV, elkold, TINY;

  if (first) {
    first = false;
    mcold = 1.0;
    PIHALF = atan(1.0) * 2.0;
    PIINV = 0.5 / PIHALF;
    elkold = PIHALF;
    TINY = 1.e-99;
  }
  m = 1.0 - mc;
  if (abs(m) < 1.d - 16) {
    elk = PIHALF;
  } else if (abs(mc - mcold) < 1.11d - 16 * mc) {
    elk = elkold;
  } else if (mc < TINY) {
    elk = 1.3862943611198906 - 0.5 * log(TINY);
  } else if (mc < 1.11d - 16) {
    elk = 1.3862943611198906 - 0.5 * log(mc);
  } else if (mc < 0.1) {
    nome=mc * (D1 + mc * (D2 + mc * (D3 + mc * (D4 + mc * (D5 + mc * (D6
          + mc * (D7 + mc * (D8 + mc * (D9 + mc * (D10 + mc * (D11 + mc * (D12
          + mc * (D13 + mc * D14)))))))))))));
    mx = mc - 0.05;
    //
    //  K'
    //
    kkc = 1.591003453790792180 + mx * (
         0.416000743991786912 + mx * (
         0.245791514264103415 + mx * (
         0.179481482914906162 + mx * (
         0.144556057087555150 + mx * (
         0.123200993312427711 + mx * (
         0.108938811574293531 + mx * (
         0.098853409871592910 + mx * (
         0.091439629201749751 + mx * (
         0.085842591595413900 + mx * (
         0.081541118718303215))))))))));
    elk = -kkc * PIINV * log(nome)
  } else if (m <= 0.1) {
    mx = m - 0.05;
    elk = 1.591003453790792180 + mx * (
         0.416000743991786912 + mx * (
         0.245791514264103415 + mx * (
         0.179481482914906162 + mx * (
         0.144556057087555150 + mx * (
         0.123200993312427711 + mx * (
         0.108938811574293531 + mx * (
         0.098853409871592910 + mx * (
         0.091439629201749751 + mx * (
         0.085842591595413900 + mx * (
         0.081541118718303215))))))))));
  } else if (m <= 0.2) {
    mx = m - 0.15;
    elk = 1.635256732264579992 + mx * (
         0.471190626148732291 + mx * (
         0.309728410831499587 + mx * (
         0.252208311773135699 + mx * (
         0.226725623219684650 + mx * (
         0.215774446729585976 + mx * (
         0.213108771877348910 + mx * (
         0.216029124605188282 + mx * (
         0.223255831633057896 + mx * (
         0.234180501294209925 + mx * (
         0.248557682972264071 + mx * (
         0.266363809892617521 + mx * (
         0.287728452156114668))))))))))));
  } else if (m <= 0.3) {
    mx = m - 0.25;
    elk = 1.685750354812596043 + mx * (
         0.541731848613280329 + mx * (
         0.401524438390690257 + mx * (
         0.369642473420889090 + mx * (
         0.376060715354583645 + mx * (
         0.405235887085125919 + mx * (
         0.453294381753999079 + mx * (
         0.520518947651184205 + mx * (
         0.609426039204995055 + mx * (
         0.724263522282908870 + mx * (
         0.871013847709812357 + mx * (
         1.057652872753547036)))))))))));
  } else if (m <= 0.4) {
    mx = m - 0.35;
    elk = 1.744350597225613243 + mx * (
         0.634864275371935304 + mx * (
         0.539842564164445538 + mx * (
         0.571892705193787391 + mx * (
         0.670295136265406100 + mx * (
         0.832586590010977199 + mx * (
         1.073857448247933265 + mx * (
         1.422091460675497751 + mx * (
         1.920387183402304829 + mx * (
         2.632552548331654201 + mx * (
         3.652109747319039160 + mx * (
         5.115867135558865806 + mx * (
         7.224080007363877411))))))))))));
  } else if (m <= 0.5) {
    mx = m - 0.45;
    elk = 1.813883936816982644 + mx * (
         0.763163245700557246 + mx * (
         0.761928605321595831 + mx * (
         0.951074653668427927 + mx * (
         1.315180671703161215 + mx * (
         1.928560693477410941 + mx * (
         2.937509342531378755 + mx * (
         4.594894405442878062 + mx * (
         7.330071221881720772 + mx * (
         11.87151259742530180 + mx * (
         19.45851374822937738 + mx * (
         32.20638657246426863 + mx * (
         53.73749198700554656 + mx * (
         90.27388602940998849)))))))))))));
  } else if (m <= 0.6) {
    mx = m - 0.55;
    elk = 1.898924910271553526 + mx * (
         0.950521794618244435 + mx * (
         1.151077589959015808 + mx * (
         1.750239106986300540 + mx * (
         2.952676812636875180 + mx * (
         5.285800396121450889 + mx * (
         9.832485716659979747 + mx * (
         18.78714868327559562 + mx * (
         36.61468615273698145 + mx * (
         72.45292395127771801 + mx * (
         145.1079577347069102 + mx * (
         293.4786396308497026 + mx * (
         598.3851815055010179 + mx * (
         1228.420013075863451 + mx * (
         2536.529755382764488))))))))))))));
  } else if (m <= 0.7) {
    mx = m - 0.65;
    elk = 2.007598398424376302 + mx * (
         1.248457231212347337 + mx * (
         1.926234657076479729 + mx * (
         3.751289640087587680 + mx * (
         8.119944554932045802 + mx * (
         18.66572130873555361 + mx * (
         44.60392484291437063 + mx * (
         109.5092054309498377 + mx * (
         274.2779548232413480 + mx * (
         697.5598008606326163 + mx * (
         1795.716014500247129 + mx * (
         4668.381716790389910 + mx * (
         12235.76246813664335 + mx * (
         32290.17809718320818 + mx * (
         85713.07608195964685 + mx * (
         228672.1890493117096 + mx * (
         612757.2711915852774))))))))))))))));
  } else if (m <= 0.8) {
    mx = m - 0.75;
    elk = 2.156515647499643235 + mx * (
         1.791805641849463243 + mx * (
         3.826751287465713147 + mx * (
         10.38672468363797208 + mx * (
         31.40331405468070290 + mx * (
         100.9237039498695416 + mx * (
         337.3268282632272897 + mx * (
         1158.707930567827917 + mx * (
         4060.990742193632092 + mx * (
         14454.00184034344795 + mx * (
         52076.66107599404803 + mx * (
         189493.6591462156887 + mx * (
         695184.5762413896145 + mx * (
         2.567994048255284686d6 + mx * (
         9.541921966748386322d6 + mx * (
         3.563492744218076174d7 + mx * (
         1.336692984612040871d8 + mx * (
         5.033521866866284541d8 + mx * (
         1.901975729538660119d9 + mx * (
         7.208915015330103756d9)))))))))))))))))));
  } else if (m <= 0.85) {
    mx = m - 0.825;
    elk = 2.318122621712510589 + mx * (
         2.616920150291232841 + mx * (
         7.897935075731355823 + mx * (
         30.50239715446672327 + mx * (
         131.4869365523528456 + mx * (
         602.9847637356491617 + mx * (
         2877.024617809972641 + mx * (
         14110.51991915180325 + mx * (
         70621.44088156540229 + mx * (
         358977.2665825309926 + mx * (
         1.847238263723971684d6 + mx * (
         9.600515416049214109d6 + mx * (
         5.030767708502366879d7 + mx * (
         2.654441886527127967d8 + mx * (
         1.408862325028702687d9 + mx * (
         7.515687935373774627d9)))))))))))))));
  } else {
    mx = m - 0.875;
    elk = 2.473596173751343912 + mx * (
         3.727624244118099310 + mx * (
         15.60739303554930496 + mx * (
         84.12850842805887747 + mx * (
         506.9818197040613935 + mx * (
         3252.277058145123644 + mx * (
         21713.24241957434256 + mx * (
         149037.0451890932766 + mx * (
         1.043999331089990839d6 + mx * (
         7.427974817042038995d6 + mx * (
         5.350383967558661151d7 + mx * (
         3.892498869948708474d8 + mx * (
         2.855288351100810619d9 + mx * (
         2.109007703876684053d10 + mx * (
         1.566998339477902014d11 + mx * (
         1.170222242422439893d12 + mx * (
         8.777948323668937971d12 + mx * (
         6.610124275248495041d13 + mx * (
         4.994880537133887989d14 + mx * (
         3.785974339724029920d15)))))))))))))))))));
  }

  mcold = mc;
  elkold = elk;
}

}  // namespace numerics
}  // namespace principia
