(* ::Package:: *)

IntegrationErrorPlot[errorData_,names_]:=
With[{
  series = Length[names],
  errorDataAndEmpty = Append[errorData,{0,0}],
  minWork = Min[First/@Flatten[errorData, 1]],
  maxWork = Max[First/@Flatten[errorData, 1]],
  colour = ColorData[97]},
  DynamicModule[
    {plots =
       Function[
         i,
         ListLogLogPlot[
           errorDataAndEmpty[[i]],
           PlotRange->{{minWork, maxWork}, {1*^-17, 1}},
           ImageSize->1200,
           AxesOrigin->{50, 2^-53},
           PlotStyle->{colour[i], PointSize[0.001]}]]/@Range[series+1],
     visible = Append[ConstantArray[False, series], True]},
    Dynamic[
      Row[
      {Show[plots[[Select[Range[series+1], visible[[#]]&]]]],
       SwatchLegend[
         If[visible[[#]], colour[#], Transparent]&/@Range[series],
         Toggler[Dynamic[visible[[#]]], {True->names[[#]], False->names[[#]]}]&/@(Range[series])]}]],
    SaveDefinitions->True]]
