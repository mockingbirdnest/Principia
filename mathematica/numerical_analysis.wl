(* ::Package:: *)

IntegrationErrorPlot[
  errorData_,
  names_,
  errorKind_,
  upperError_: 1,
  lowerRelative_: 1.*^-17] :=
 With[
  {series = Length[names],
   errorDataAndEmpty = Append[errorData, {0, 0}],
   minWork = Min[First /@ Flatten[errorData, 1]],
   maxWork = Max[First /@ Flatten[errorData, 1]],
   unit = QuantityUnit@UnitSimplify@Last@Last@Last@errorData,
   (*The Mathematica 10 default indexed colour scheme*)
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
       PrintTemporary["Rendering " <> errorKind <> ":"];
       (*Beware the Jabberwock, my son! |temporaryCell| must be deleted
         before we leave the DynamicModule, otherwise it will refer to the
         out of scope variable |currentTasks|.*)
       temporaryCell = PrintTemporary[Dynamic[Column[currentTasks]]];
       Function[
         i,
         If[i <= series,
          currentTasks[[$KernelID /. kernels]] =
           names[[i]] <> " (" <> ToString[i] <> "/" <>
            ToString[series] <> ")"];
         ListLogLogPlot[
          errorDataAndEmpty[[i]],
          PlotRange -> {{minWork, maxWork},
                         {lowerRelative upperError, upperError}},
          ImageSize -> 1200,
          AxesLabel -> {"Evaluations", errorKind <>
             " (" <> ToString[unit, TraditionalForm] <> ")"},
          PlotStyle -> {colour[i], PointSize[0.001]}]]~ParallelMap~
        Range[series + 1]),
     visible = Append[ConstantArray[False, series], True]},
    NotebookDelete@temporaryCell;
    PrintTemporary["Done."];
    Dynamic[
     Row[
      {Show[plots[[Select[Range[series + 1], visible[[#]] &]]]],
       SwatchLegend[
        If[visible[[#]], colour[#], Transparent] & /@ Range[series],
        Toggler[
           Dynamic[visible[[#]]],
           {True -> names[[#]], False -> names[[#]]}] & /@
         Range[series]]}]],
    SaveDefinitions -> True]]]
