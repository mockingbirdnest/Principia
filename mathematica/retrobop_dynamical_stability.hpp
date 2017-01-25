#pragma once

namespace principia {
namespace mathematica {

// Generates data for plotting the positions of the Joolian system over 5 years.
// See |AnalyseGlobalError| for a motivation for that figure.
void PlotPredictableYears();

// Generates data for apsis plots of the Joolian moons over a century, as well
// as various properties of Bop's orbit, and Bop's separations to Pol and Tylo.
void PlotCentury();

// Uses an integration of the system at the base |step|, and integration of the
// system at |step / 2|, and a slew of integrations of millimetrically perturbed
// systems at |step| to show that the error from numerical integration
// (estimated in forward error as the error between the integration at |step|
// and the one at |step / 2|) corresponds to a submillimetric backward error
// (as the millimetrically perturbed cluster diverges much faster than the
// refined integration does).
// Moreover, shows that the forward error from millimetric perturbations is less
// than 1e6 m after 5 a: this duration is deemed "predictable" (see above).
void AnalyseGlobalError();

// Integrates a cluster of millimetrically perturbed systems, discarding them as
// the integrations encounter high local errors, currently very badly estimated,
// or eject Joolian moons.
// TODO(egg): This should be reworked to actually estimate the local error and
// potentially control it, and to detect close encounters.
void StatisticallyAnalyseStability();

}  // namespace mathematica
}  // namespace principia
