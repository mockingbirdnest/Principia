(* ::Package:: *)

IntegrationErrorPlot[errorData_, names_, errorKind_] :=
With[
  {series = Length[names],
   errorDataAndEmpty = Append[errorData,{0,0}],
   minWork = Min[First/@Flatten[errorData, 1]],
   maxWork = Max[First/@Flatten[errorData, 1]],
   unit = QuantityUnit@UnitSimplify@Last@Last@Last@errorData,
   colour = ColorData[97]},
  Module[
    {kernels = ParallelTable[$KernelID -> i, {i, $KernelCount}],
     currentTasks,
     temporaryCell},
    DynamicModule[
      {plots = (
           SetSharedVariable[kernels];
           SetSharedVariable[currentTasks];
           currentTasks = ConstantArray["", Length[kernels]];
           PrintTemporary["Rendering "<>errorKind<>":"];
           temporaryCell = PrintTemporary[Dynamic[Column[currentTasks]]];
           Function[
             i,
	         If[i <= series,
                currentTasks[[$KernelID /. kernels]] =
                  names[[i]]<>" ("<>ToString[i]<>"/"<>ToString[series]<>")"];
             ListLogLogPlot[
               errorDataAndEmpty[[i]],
               PlotRange -> {{minWork, maxWork}, {1*^-17, 1}},
               ImageSize -> 1200,
               AxesLabel -> {"Evaluations", errorKind<>
                               " ("<>ToString[unit, TraditionalForm]<>")"},
               PlotStyle -> {colour[i], PointSize[0.001]}]]~ParallelMap~Range[series+1]),
       visible = Append[ConstantArray[False, series], True]},
      NotebookDelete@temporaryCell;
      PrintTemporary["Done."]
      Dynamic[
        Row[
        {Show[plots[[Select[Range[series+1], visible[[#]]&]]]],
         SwatchLegend[
           If[visible[[#]], colour[#], Transparent]&/@Range[series],
           Toggler[Dynamic[visible[[#]]],
                   {True -> names[[#]], False -> names[[#]]}]&/@(Range[series])]}]],
      SaveDefinitions->True]]]
