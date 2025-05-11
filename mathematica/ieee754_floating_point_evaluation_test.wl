(* ::Package:: *)

(* ::Input::Initialization:: *)
Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point_evaluation.wl"}]]


(* ::Input::Initialization:: *)
On[Assert]


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary64]
SetRoundingMode[NearestTiesToEven]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1]==1]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1.5]==3/2]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[0.1`100]==CorrectlyRound[1/10]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3]==CorrectlyRound[1/3]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1+2]==3]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3+1/3]==2CorrectlyRound[1/3]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3+1/5]==CorrectlyRound[CorrectlyRound[1/3]+CorrectlyRound[1/5]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3-1/5]==CorrectlyRound[CorrectlyRound[1/3]-CorrectlyRound[1/5]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)*(1/5)]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)^2]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)^3]==
CorrectlyRound[
CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[1/3]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)^4]==
CorrectlyRound[
CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[N[1/3,100]]==CorrectlyRound[1/3]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[\[Pi]]=!=CorrectlyRound[\[Pi]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]+CorrectlyRound[1/7]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]==CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]+CorrectlyRound[1/7]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]=!=IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]]
