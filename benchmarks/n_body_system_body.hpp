#pragma once

#include "geometry/grassmann.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace benchmarks {

void SolarSystem1643To1727() {
  using physics::Body;
  using physics::NBodySystem;
  using quantities::GravitationalParameter;
  using quantities::Pow;
  using quantities::SIUnit;
  using quantities::Time;
  using si::Day;
  using si::Kilo;
  using si::Kilogram;
  using si::Metre;
  using si::Second;

  typedef Vector<Length, ICRFJ2000EclipticFrame> Position;
  typedef Vector<Speed, ICRFJ2000EclipticFrame> Velocity;

  // Number of days since the JD epoch. JD2436116.3115 is the time of the launch
  // of Простейший Спутник-1.
  Time t0 = 2436116.3115 * Day;

  // Star.
  Body<ICRFJ2000EclipticFrame> sun(
      1.3271244004193938E+11 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  sun.AppendToTrajectory(
      Position({ 1.138350928138014E+06 * Kilo(Metre),
                 6.177753685036716E+05 * Kilo(Metre),
                -3.770941657504326E+04 * Kilo(Metre)}),
      Velocity({-5.067456621846211E-03 * Kilo(Metre) / Second,
                 1.259599196445122E-02 * Kilo(Metre) / Second,
                 9.778588606052481E-05 * Kilo(Metre) / Second}),
      t0);

  // Planets.

  // Gas giants.
  Body<ICRFJ2000EclipticFrame> jupiter(
      126686511 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  jupiter.AppendToTrajectory(
      Position({-7.950209667306620E+08 * Kilo(Metre),
                -1.784285526424396E+08 * Kilo(Metre),
                 1.853825132237791E+07 * Kilo(Metre)}),
      Velocity({ 2.709330231918198E+00 * Kilo(Metre) / Second,
                -1.213073724288562E+01 * Kilo(Metre) / Second,
                -1.088748435062713E-02 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> saturn(
      37931207.8 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  saturn.AppendToTrajectory(
      Position({-3.774715321901159E+08 * Kilo(Metre),
                -1.451892263379818E+09 * Kilo(Metre),
                 4.040621083792380E+07 * Kilo(Metre)}),
      Velocity({ 8.817029873536633E+00 * Kilo(Metre) / Second,
                -2.466058486223613E+00 * Kilo(Metre) / Second,
                -3.068419809533604E-01 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> neptune(
      6835107 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  neptune.AppendToTrajectory(
      Position({-3.810689792831146E+09 * Kilo(Metre),
                -2.456423858579051E+09 * Kilo(Metre),
                 1.383694320077938E+08 * Kilo(Metre)}),
      Velocity({ 2.913267720085410E+00 * Kilo(Metre) / Second,
                -4.535247383721019E+00 * Kilo(Metre) / Second,
                 2.589759251085161E-02 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> uranus(
      5793966 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  uranus.AppendToTrajectory(
      Position({-1.729995609344851E+09 * Kilo(Metre),
                 2.159967050539728E+09 * Kilo(Metre),
                 3.048735047038063E+07 * Kilo(Metre)}),
      Velocity({-5.366539669972795E+00 * Kilo(Metre) / Second,
                -4.575802196749351E+00 * Kilo(Metre) / Second,
                 5.261322980347850E-02 * Kilo(Metre) / Second}),
      t0);

  // Telluric planets.
  Body<ICRFJ2000EclipticFrame> earth(
      398600.440 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  earth.AppendToTrajectory(
      Position({ 1.475150112055673E+08 * Kilo(Metre),
                 3.144435102288270E+07 * Kilo(Metre),
                -3.391764309344300E+04 * Kilo(Metre)}),
      Velocity({-6.635753510543799E+00 * Kilo(Metre) / Second,
                 2.904321639216012E+01 * Kilo(Metre) / Second,
                 3.125252418990812E-03 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> venus(
      324858.63 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  venus.AppendToTrajectory(
      Position({ 6.084974577091119E+07 * Kilo(Metre),
                -9.037413730207849E+07 * Kilo(Metre),
                -4.719158908401959E+06 * Kilo(Metre)}),
      Velocity({ 2.903958257174759E+01 * Kilo(Metre) / Second,
                 1.910383147602264E+01 * Kilo(Metre) / Second,
                -1.418780340302349E+00 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> mars(
      42828.3 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  mars.AppendToTrajectory(
      Position({-2.440047184660406E+08 * Kilo(Metre),
                -2.002994580992744E+07 * Kilo(Metre),
                 5.577600092368793E+06 * Kilo(Metre)}),
      Velocity({ 2.940381268511949E+00 * Kilo(Metre) / Second,
                -2.206625841382794E+01 * Kilo(Metre) / Second,
                -5.348179460834037E-01 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame> mercury(
      22032.09 * Pow<3>(Kilo(Metre)) / Pow<2>(Second));
  mercury.AppendToTrajectory(
      Position({-3.013851560892715E+07 * Kilo(Metre),
                 3.823388939456400E+07 * Kilo(Metre),
                 5.907240907643730E+06 * Kilo(Metre)}),
      Velocity({-4.731017449071709E+01 * Kilo(Metre) / Second,
                -2.918747853895398E+01 * Kilo(Metre) / Second,
                 1.963450229872517E+00 * Kilo(Metre) / Second}),
      t0);

  // End of planets.


}

}
}
