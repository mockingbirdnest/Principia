#pragma once

#include "testing_utilities/solar_system.hpp"

#include <vector>

#include "geometry/grassmann.hpp"
#include "physics/body.hpp"
#include "physics/n_body_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace testing_utilities {

physics::NBodySystem<ICRFJ2000EclipticFrame> SolarSystemAtSputnikLaunch() {
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

  // All data is from the Jet Propulsion Laboratory's HORIZONS system unless
  // otherwise specified.

  // Star.
  Body<ICRFJ2000EclipticFrame>* sun(
      new Body<ICRFJ2000EclipticFrame>(
          1.3271244004193938E+11 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  sun->AppendToTrajectory(
      Position({ 1.138350928138014E+06 * Kilo(Metre),
                 6.177753685036716E+05 * Kilo(Metre),
                -3.770941657504326E+04 * Kilo(Metre)}),
      Velocity({-5.067456621846211E-03 * Kilo(Metre) / Second,
                 1.259599196445122E-02 * Kilo(Metre) / Second,
                 9.778588606052481E-05 * Kilo(Metre) / Second}),
      t0);

  // Planets.

  // Gas giants.
  Body<ICRFJ2000EclipticFrame>* jupiter(
      new Body<ICRFJ2000EclipticFrame>(
          126686511 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  jupiter->AppendToTrajectory(
      Position({-7.950209667306620E+08 * Kilo(Metre),
                -1.784285526424396E+08 * Kilo(Metre),
                 1.853825132237791E+07 * Kilo(Metre)}),
      Velocity({ 2.709330231918198E+00 * Kilo(Metre) / Second,
                -1.213073724288562E+01 * Kilo(Metre) / Second,
                -1.088748435062713E-02 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* saturn(
      new Body<ICRFJ2000EclipticFrame>(
          37931207.8 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  saturn->AppendToTrajectory(
      Position({-3.774715321901159E+08 * Kilo(Metre),
                -1.451892263379818E+09 * Kilo(Metre),
                 4.040621083792380E+07 * Kilo(Metre)}),
      Velocity({ 8.817029873536633E+00 * Kilo(Metre) / Second,
                -2.466058486223613E+00 * Kilo(Metre) / Second,
                -3.068419809533604E-01 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* neptune(
      new Body<ICRFJ2000EclipticFrame>(
          6835107 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  neptune->AppendToTrajectory(
      Position({-3.810689792831146E+09 * Kilo(Metre),
                -2.456423858579051E+09 * Kilo(Metre),
                 1.383694320077938E+08 * Kilo(Metre)}),
      Velocity({ 2.913267720085410E+00 * Kilo(Metre) / Second,
                -4.535247383721019E+00 * Kilo(Metre) / Second,
                 2.589759251085161E-02 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* uranus(
      new Body<ICRFJ2000EclipticFrame>(
          5793966 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  uranus->AppendToTrajectory(
      Position({-1.729995609344851E+09 * Kilo(Metre),
                 2.159967050539728E+09 * Kilo(Metre),
                 3.048735047038063E+07 * Kilo(Metre)}),
      Velocity({-5.366539669972795E+00 * Kilo(Metre) / Second,
                -4.575802196749351E+00 * Kilo(Metre) / Second,
                 5.261322980347850E-02 * Kilo(Metre) / Second}),
      t0);

  // Telluric planets.
  Body<ICRFJ2000EclipticFrame>* earth(
      new Body<ICRFJ2000EclipticFrame>(
          398600.440 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  earth->AppendToTrajectory(
      Position({ 1.475150112055673E+08 * Kilo(Metre),
                 3.144435102288270E+07 * Kilo(Metre),
                -3.391764309344300E+04 * Kilo(Metre)}),
      Velocity({-6.635753510543799E+00 * Kilo(Metre) / Second,
                 2.904321639216012E+01 * Kilo(Metre) / Second,
                 3.125252418990812E-03 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* venus(
      new Body<ICRFJ2000EclipticFrame>(
          324858.63 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  venus->AppendToTrajectory(
      Position({ 6.084974577091119E+07 * Kilo(Metre),
                -9.037413730207849E+07 * Kilo(Metre),
                -4.719158908401959E+06 * Kilo(Metre)}),
      Velocity({ 2.903958257174759E+01 * Kilo(Metre) / Second,
                 1.910383147602264E+01 * Kilo(Metre) / Second,
                -1.418780340302349E+00 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* mars(
      new Body<ICRFJ2000EclipticFrame>(
          42828.3 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  mars->AppendToTrajectory(
      Position({-2.440047184660406E+08 * Kilo(Metre),
                -2.002994580992744E+07 * Kilo(Metre),
                 5.577600092368793E+06 * Kilo(Metre)}),
      Velocity({ 2.940381268511949E+00 * Kilo(Metre) / Second,
                -2.206625841382794E+01 * Kilo(Metre) / Second,
                -5.348179460834037E-01 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* mercury(
      new Body<ICRFJ2000EclipticFrame>(
          22032.09 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  mercury->AppendToTrajectory(
      Position({-3.013851560892715E+07 * Kilo(Metre),
                 3.823388939456400E+07 * Kilo(Metre),
                 5.907240907643730E+06 * Kilo(Metre)}),
      Velocity({-4.731017449071709E+01 * Kilo(Metre) / Second,
                -2.918747853895398E+01 * Kilo(Metre) / Second,
                 1.963450229872517E+00 * Kilo(Metre) / Second}),
      t0);

  // End of planets.

  // Satellite of Jupiter.
  Body<ICRFJ2000EclipticFrame>* ganymede(
      new Body<ICRFJ2000EclipticFrame>(1482E20 * Kilogram));
  ganymede->AppendToTrajectory(
      Position({-7.942681422941415E+08 * Kilo(Metre),
                -1.776681035234876E+08 * Kilo(Metre),
                 1.857215495334835E+07 * Kilo(Metre)}),
      Velocity({-5.026319376504355E+00 * Kilo(Metre) / Second,
                -4.481735740234995E+00 * Kilo(Metre) / Second,
                 1.326192167761359E-01 * Kilo(Metre) / Second}),
      t0);

  // Satellite of Saturn.
  Body<ICRFJ2000EclipticFrame>* titan(
      new Body<ICRFJ2000EclipticFrame>(
          8978.13 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  titan->AppendToTrajectory(
      Position({-3.771930512714775E+08 * Kilo(Metre),
                -1.452931696594699E+09 * Kilo(Metre),
                 4.091643033375849E+07 * Kilo(Metre)}),
      Velocity({ 1.433381483669744E+01 * Kilo(Metre) / Second,
                -1.422590492527597E+00 * Kilo(Metre) / Second,
                -1.375826555026097E+00 * Kilo(Metre) / Second}),
      t0);

  // Satellites of Jupiter.
  Body<ICRFJ2000EclipticFrame>* callisto(
      new Body<ICRFJ2000EclipticFrame>(1076E20 * Kilogram));
  callisto->AppendToTrajectory(
      Position({-7.951805452047400E+08 * Kilo(Metre),
                -1.802957437059298E+08 * Kilo(Metre),
                 1.847154088070625E+07 * Kilo(Metre)}),
      Velocity({ 1.091928199422218E+01 * Kilo(Metre) / Second,
                -1.278098875182818E+01 * Kilo(Metre) / Second,
                 5.878649120351949E-02 * Kilo(Metre) / Second}),
      t0);
  Body<ICRFJ2000EclipticFrame>* io(
      new Body<ICRFJ2000EclipticFrame>(893.3E20 * Kilogram));
  io->AppendToTrajectory(
      Position({-7.946073188298367E+08 * Kilo(Metre),
                -1.783491436977172E+08 * Kilo(Metre),
                 1.854699192614355E+07 * Kilo(Metre)}),
      Velocity({-5.049684272040893E-01 * Kilo(Metre) / Second,
                 4.916473261567652E+00 * Kilo(Metre) / Second,
                 5.469177855959977E-01 * Kilo(Metre) / Second}),
      t0);

  // Satellite of Earth.
  Body<ICRFJ2000EclipticFrame>* moon(
      new Body<ICRFJ2000EclipticFrame>(
          4902.798 * Pow<3>(Kilo(Metre)) / Pow<2>(Second)));
  moon->AppendToTrajectory(
      Position({ 1.478545271460863E+08 * Kilo(Metre),
                 3.122566749814625E+07 * Kilo(Metre),
                 1.500491219719345E+03 * Kilo(Metre)}),
      Velocity({-6.099833968412930E+00 * Kilo(Metre) / Second,
                 2.985006033154299E+01 * Kilo(Metre) / Second,
                -1.952438319420470E-02 * Kilo(Metre) / Second}),
      t0);

  // Satellite of Jupiter.
  Body<ICRFJ2000EclipticFrame>* europa(
      new Body<ICRFJ2000EclipticFrame>(479.7E20 * Kilogram));
  europa->AppendToTrajectory(
      Position({-7.944180333947762E+08 * Kilo(Metre),
                -1.787346439588362E+08 * Kilo(Metre),
                 1.853675837527557E+07 * Kilo(Metre)}),
      Velocity({ 8.811255547505889E+00 * Kilo(Metre) / Second,
                 5.018147960240774E-02 * Kilo(Metre) / Second,
                 6.162195631257494E-01 * Kilo(Metre) / Second}),
      t0);

  // Satellite of Neptune.
  Body<ICRFJ2000EclipticFrame>* triton(
      new Body<ICRFJ2000EclipticFrame>(214.7E20 * Kilogram));
  triton->AppendToTrajectory(
      Position({-3.810797098554279E+09 * Kilo(Metre),
                -2.456691608348630E+09 * Kilo(Metre),
                 1.381629136719314E+08 * Kilo(Metre)}),
      Velocity({-1.047462448797063E+00 * Kilo(Metre) / Second,
                -4.404556713303486E+00 * Kilo(Metre) / Second,
                 1.914469843538767E+00 * Kilo(Metre) / Second}),
      t0);

  // Dwarf planet (scattered disc object).
  // Mass from Brown, Michael E.; Schaller, Emily L. (15 June 2007).
  // "The Mass of Dwarf Planet Eris", in Science, through Wikipedia.
  Body<ICRFJ2000EclipticFrame>* eris(
      new Body<ICRFJ2000EclipticFrame>(1.67E22 * Kilogram));
  eris->AppendToTrajectory(
      Position({ 1.317390066862979E+10 * Kilo(Metre),
                 2.221403321600002E+09 * Kilo(Metre),
                -5.736076877456254E+09 * Kilo(Metre)}),
      Velocity({ 4.161883594267296E-01 * Kilo(Metre) / Second,
                 1.872714752602233E+00 * Kilo(Metre) / Second,
                 1.227093842948539E+00 * Kilo(Metre) / Second}),
      t0);

  // Dwarf planet (Kuiper belt object).
  Body<ICRFJ2000EclipticFrame>* pluto(
      new Body<ICRFJ2000EclipticFrame>(1.307E22 * Kilogram));
  pluto->AppendToTrajectory(
      Position({-4.406985590968750E+09 * Kilo(Metre),
                 2.448731153209013E+09 * Kilo(Metre),
                 1.012525975599311E+09 * Kilo(Metre)}),
      Velocity({-1.319871918266467E+00 * Kilo(Metre) / Second,
                -5.172112237151897E+00 * Kilo(Metre) / Second,
                 9.407707128142039E-01 * Kilo(Metre) / Second}),
      t0);

  // End of celestial bodies.

  std::vector<Body<ICRFJ2000EclipticFrame>*> const* bodies(
      new std::vector<Body<ICRFJ2000EclipticFrame>*> const{
          sun, jupiter, saturn, neptune, uranus, earth, venus, mars, mercury,
          ganymede, titan, callisto, io, moon, europa, triton, eris, pluto});
  NBodySystem<ICRFJ2000EclipticFrame> result(bodies);
}

}  // namespace testing_utilities
}  // namespace principia
