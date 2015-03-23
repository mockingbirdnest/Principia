(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
<<"numerical_analysis.wl";
Run[FileNameJoin[{"..", "Release", "mathematica.exe"}]];
<<"simple_harmonic_motion_graphs.generated.wl";
Export[
 "shm_energy_error.cdf",
 IntegrationErrorPlot[eErrorData, names, "maximal energy error"]];
Export[
 "shm_position_error.cdf",
 IntegrationErrorPlot[qErrorData, names, "maximal position error"]];
Export[
 "shm_velocity_error.cdf",
 IntegrationErrorPlot[vErrorData, names, "maximal velocity error"]];
<<"kepler_problem_graphs.generated.wl";
Export[
 "kepler_energy_error.cdf",
 IntegrationErrorPlot[eErrorData, names, "maximal energy error", 2.*^9]];
Export[
 "kepler_position_error.cdf",
 IntegrationErrorPlot[qErrorData, names, "maximal position error"]];
Export[
 "kepler_velocity_error.cdf",
 IntegrationErrorPlot[vErrorData, names, "maximal velocity error"]];
