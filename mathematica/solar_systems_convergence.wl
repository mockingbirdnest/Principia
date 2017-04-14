(* ::Package:: *)

(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"ksp_system_convergence.generated.wl"}]]*)


(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"solar_system_convergence.generated.wl"}]]*)


(* ::Input:: *)
(*getTime[ppaConvergence_]:=Map[{#[[4]],#[[2]]}&,ppaConvergence]*)


(* ::Input:: *)
(*getStep[ppaConvergence_]:=Map[{#[[1]],#[[2]]}&,ppaConvergence]*)


(* ::Input:: *)
(*blanesSolarSystemConvergence=ppaSolarSystemConvergence0*)


(* ::Input:: *)
(*mclachlanSolarSystemConvergence=ppaSolarSystemConvergence1*)


(* ::Input:: *)
(*quinlan8SolarSystemConvergence=ppaSolarSystemConvergence2*)


(* ::Input:: *)
(*quinlan10SolarSystemConvergence=ppaSolarSystemConvergence3*)


(* ::Input:: *)
(*quinlan12SolarSystemConvergence=ppaSolarSystemConvergence4*)


(* ::Input:: *)
(*preferredSolarSystemConvergence=ppaSolarSystemConvergence5*)


(* ::Input:: *)
(*Show[*)
(*ListLogLogPlot[getTime[blanesSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],*)
(*ListLogLogPlot[getTime[mclachlanSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],*)
(*ListLogLogPlot[getTime[quinlan8SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],*)
(*ListLogLogPlot[getTime[quinlan10SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],*)
(*ListLogLogPlot[getTime[quinlan12SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],*)
(*ListLogLogPlot[getTime[preferredSolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]*)
(*]*)


(* ::Input:: *)
(*Show[*)
(*ListLogLogPlot[getStep[blanesSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],*)
(*ListLogLogPlot[getStep[mclachlanSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],*)
(*ListLogLogPlot[getStep[quinlan8SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],*)
(*ListLogLogPlot[getStep[quinlan10SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],*)
(*ListLogLogPlot[getStep[quinlan12SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],*)
(*ListLogLogPlot[getStep[preferredSolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]*)
(*]*)


(* ::Input:: *)
(*blanesKSPSystemConvergence=ppaKSPSystemConvergence0*)


(* ::Input:: *)
(*mclachlanKSPSystemConvergence=ppaKSPSystemConvergence1*)


(* ::Input:: *)
(*quinlan8KSPSystemConvergence=ppaKSPSystemConvergence2*)


(* ::Input:: *)
(*quinlan10KSPSystemConvergence=ppaKSPSystemConvergence3*)


(* ::Input:: *)
(*quinlan12KSPSystemConvergence=ppaKSPSystemConvergence4*)


(* ::Input:: *)
(*preferredKSPSystemConvergence=ppaKSPSystemConvergence5*)


(* ::Input:: *)
(*Show[*)
(*ListLogLogPlot[getTime[blanesKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],*)
(*ListLogLogPlot[getTime[mclachlanKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],*)
(*ListLogLogPlot[getTime[quinlan8KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],*)
(*ListLogLogPlot[getTime[quinlan10KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],*)
(*ListLogLogPlot[getTime[quinlan12KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],*)
(*ListLogLogPlot[getTime[preferredKSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]*)
(*]*)


(* ::Input:: *)
(*Show[*)
(*ListLogLogPlot[getStep[blanesKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],*)
(*ListLogLogPlot[getStep[mclachlanKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],*)
(*ListLogLogPlot[getStep[quinlan8KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],*)
(*ListLogLogPlot[getStep[quinlan10KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],*)
(*ListLogLogPlot[getStep[quinlan12KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],*)
(*ListLogLogPlot[getStep[preferredKSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]*)
(*]*)
