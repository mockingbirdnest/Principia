(* ::Package:: *)

Get[FileNameJoin[{NotebookDirectory[],"ksp_system_convergence.generated.wl"}]]


Get[FileNameJoin[{NotebookDirectory[],"solar_system_convergence.generated.wl"}]]


getTime[ppaConvergence_]:=Map[{#[[4]],#[[2]]}&,ppaConvergence]


getStep[ppaConvergence_]:=Map[{#[[1]],#[[2]]}&,ppaConvergence]


blanesSolarSystemConvergence=ppaSolarSystemConvergence0


mclachlanSolarSystemConvergence=ppaSolarSystemConvergence1


quinlan8SolarSystemConvergence=ppaSolarSystemConvergence2


quinlan10SolarSystemConvergence=ppaSolarSystemConvergence3


quinlan12SolarSystemConvergence=ppaSolarSystemConvergence4


preferredSolarSystemConvergence=ppaSolarSystemConvergence5


Show[
ListLogLogPlot[getTime[blanesSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],
ListLogLogPlot[getTime[mclachlanSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],
ListLogLogPlot[getTime[quinlan8SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],
ListLogLogPlot[getTime[quinlan10SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],
ListLogLogPlot[getTime[quinlan12SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],
ListLogLogPlot[getTime[preferredSolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]
]


Show[
ListLogLogPlot[getStep[blanesSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],
ListLogLogPlot[getStep[mclachlanSolarSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],
ListLogLogPlot[getStep[quinlan8SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],
ListLogLogPlot[getStep[quinlan10SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],
ListLogLogPlot[getStep[quinlan12SolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],
ListLogLogPlot[getStep[preferredSolarSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]
]


blanesKSPSystemConvergence=ppaKSPSystemConvergence0


mclachlanKSPSystemConvergence=ppaKSPSystemConvergence1


quinlan8KSPSystemConvergence=ppaKSPSystemConvergence2


quinlan10KSPSystemConvergence=ppaKSPSystemConvergence3


quinlan12KSPSystemConvergence=ppaKSPSystemConvergence4


preferredKSPSystemConvergence=ppaKSPSystemConvergence5


Show[
ListLogLogPlot[getTime[blanesKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],
ListLogLogPlot[getTime[mclachlanKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],
ListLogLogPlot[getTime[quinlan8KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],
ListLogLogPlot[getTime[quinlan10KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],
ListLogLogPlot[getTime[quinlan12KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],
ListLogLogPlot[getTime[preferredKSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]
]


Show[
ListLogLogPlot[getStep[blanesKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["red"]],
ListLogLogPlot[getStep[mclachlanKSPSystemConvergence], Joined->True,PlotStyle->RGBColor["green"]],
ListLogLogPlot[getStep[quinlan8KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["blue"]],
ListLogLogPlot[getStep[quinlan10KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["cyan"]],
ListLogLogPlot[getStep[quinlan12KSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["magenta"]],
ListLogLogPlot[getStep[preferredKSPSystemConvergence], Joined->True,PlotRange->All,PlotStyle->RGBColor["black"]]
]
