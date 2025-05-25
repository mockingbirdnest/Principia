(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


(* ::Input::Initialization:: *)
Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]];Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point_evaluation.wl"}]];


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary64]
SetRoundingMode[NearestTiesToEven]


(* ::Input::Initialization:: *)
On[Assert]


(* ::Input::Initialization:: *)
?"IEEE754FloatingPointEvaluation`*"


(* ::Section::Closed:: *)
(*IEEEEvaluate*)


(* ::Input::Initialization:: *)
?UseFMA


(* ::Input::Initialization:: *)
?IEEEEvaluate


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1]==1];
Assert[IEEEEvaluate[1+2]==3];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1.5]==3/2];
Assert[IEEEEvaluate[0.1`100]==CorrectlyRound[1/10]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3]==CorrectlyRound[1/3]];
Assert[IEEEEvaluate[1/3+1/3]==2CorrectlyRound[1/3]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3+1/5]==CorrectlyRound[CorrectlyRound[1/3]+CorrectlyRound[1/5]]];
Assert[IEEEEvaluate[1/3-1/5]==CorrectlyRound[CorrectlyRound[1/3]-CorrectlyRound[1/5]]];
Assert[IEEEEvaluate[(1/3)*(1/5)]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)^2]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]];
Assert[IEEEEvaluate[(1/3)^3]==
CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[1/3]]];
Assert[IEEEEvaluate[(1/3)^4]==
CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[N[1/3,100]]==CorrectlyRound[1/3]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[\[Pi]]==$Failed]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]+CorrectlyRound[1/7]]];
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]==CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]+CorrectlyRound[1/7]]];
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]!=IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((1/10+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10]==
Nest[CorrectlyRound[#+CorrectlyRound[1/10]]&,0,10]];
Assert[
IEEEEvaluate[((((((((1/10+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10]!=1];
Assert[
IEEEEvaluate[1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10]==$Failed];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((1/10*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10]==
Nest[CorrectlyRound[#*CorrectlyRound[1/10]]&,1,10]];
Assert[
IEEEEvaluate[((((((((1/10*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10]!=10^-10];
Assert[
IEEEEvaluate[1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10]==$Failed];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((-1/10-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10]==
Nest[CorrectlyRound[#-CorrectlyRound[1/10]]&,0,10]];
Assert[
IEEEEvaluate[((((((((-1/10-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10]!=-1];
Assert[
IEEEEvaluate[-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10]==$Failed];


(* ::Section:: *)
(*IEEEEvaluateInterval*)


(* ::Input::Initialization:: *)
?IEEEEvaluateWithAbsoluteError


(* ::Input::Initialization:: *)
$ContextPath=Prepend[$ContextPath,"IEEE754FloatingPointEvaluation`Private`"]


(* ::Input::Initialization:: *)
Assert[halfULPAbove1==2^-53];


(* ::Input::Initialization:: *)
Assert[intervalMin[Interval[{-2,3}]]==Min[Interval[{-2,3}]]];
Assert[intervalMax[Interval[{-2,3}]]==Max[Interval[{-2,3}]]];
Assert[intervalMin[Undefined]=!=Min[Undefined]];
Assert[intervalMax[Undefined]=!=Max[Undefined]];


(* ::Input::Initialization:: *)
Assert[halfULP[Interval[{1.2,1.3}]]==2^-53];
Assert[halfULP[Interval[{-2.4,1.3}]]==2^-52];
Assert[halfULP[Interval[{1,1}]]==2^-53];
Assert[halfULP[Interval[{-1/2,1/2}]]==2^-54];


(* ::Input::Initialization:: *)
Assert[addHalfULPInterval[Interval[{3/2,5/3}]]==Interval[{3/2-2^-53,5/3+2^-53}]];
Assert[addHalfULPInterval[Interval[{3/2,5}]]==Interval[{3/2-2^-51,5+2^-51}]];
Assert[addHalfULPInterval[Interval[{-5,3/2}]]==Interval[{-5-2^-51,3/2+2^-51}]];


(* ::Input:: *)
(*Assert[addHalfULPInterval[nonnegative[Interval[{3/2,5/3}]]]==nonnegative[Interval[{3/2-2^-53,5/3+2^-53}]]];*)
(*Assert[addHalfULPInterval[nonnegative[Interval[{0,5}]]]==nonnegative[Interval[{0,5+2^-51}]]];*)
(*Assert[addHalfULPInterval[nonnegative[Interval[{2^-53,5}]]]==nonnegative[Interval[{0,5+2^-51}]]];*)


(* ::Input::Initialization:: *)
$ContextPath=Drop[$ContextPath,1]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[2*3+4,UseFMA->True]=={Interval[{10,10}],Interval[{-2^-50,2^-50}]}];
Assert[IEEEEvaluateWithAbsoluteError[2*3+4,UseFMA->False]=={Interval[{10,10}],Interval[{-2^-51-2^-50,2^-51+2^-50}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[1+2]=={Interval[{3,3}],Interval[{-2^-52,2^-52}]}];
Assert[IEEEEvaluateWithAbsoluteError[2*3]=={Interval[{6,6}],Interval[{-2^-51,2^-51}]}];
Assert[IEEEEvaluateWithAbsoluteError[2/3]=={Interval[{CorrectlyRound[2/3],CorrectlyRound[2/3]}],Interval[{-2^-54,2^-54}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[-Interval[{-2,1}]]=={Interval[{-1,2}],Interval[{0,0}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,2}]^2]=={Interval[{0,4}],Interval[{-2^-51,2^-51}]}];
(* The error after the squaring is [0, 2^-50]. *)
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^3]=={Interval[{-1,27}],Interval[{-3 2^-50-2^-49,3 2^-50+2^-49}]}];
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^4]=={Interval[{0,81}],Interval[{-18 2^-50-2^-47,2^-100+18 2^-50+2^-47}]}];


(* ::Input:: *)
(*IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^4]*)


(* ::Input:: *)
(*{v,\[Delta]}=IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^2]*)


(* ::Input:: *)
(*2 v \[Delta]+\[Delta]^2+Interval[{-2^-47,2^-47}]*)


(* ::Input:: *)
(*-18 2^-50-2^-47*)


(* ::Input:: *)
(*2^-50*)


(* ::Input:: *)
(*3 2^-50+2^-49*)


(* ::Input:: *)
(*NumberQ[Interval[{-1,2}]]*)


(* ::Input:: *)
(*Trace[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^3],applyOp|Expand|ReplaceAll|Interval]*)


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[1]==1];
Assert[IEEEEvaluateWithAbsoluteError[0.1]==CorrectlyRound[0.1]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[Undefined]===Undefined];
Assert[IEEEEvaluateWithAbsoluteError[Interval[{1,Undefined}]]===Interval[{1,Undefined}]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[2,3]==$Failed];


(* ::Input:: *)
(*(Interval[{CorrectlyRound[intervalMin[hh^2]],CorrectlyRound[intervalMax[hh^2]]}]-hh)/.hh->aa*)


(* ::Input:: *)
(*IEEEEvaluateWithAbsoluteError[3 Interval[{-1,2}]^2]*)


(* ::Input:: *)
(*aa=Interval[{-1,2}]*)


(* ::Input:: *)
(*IEEEEvaluateWithAbsoluteError[3 aa^2]*)


(* ::Input:: *)
(*PowerExpand[IEEEEvaluateWithAbsoluteError[3 hh^2]/.hh->aa]*)


(* ::Input:: *)
(*Trace[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^3],IEEE754FloatingPointEvaluation`Private`addHalfULPInterval|IEEE754FloatingPointEvaluation`Private`evae|IEEE754FloatingPointEvaluation`Private`halfULP]*)


(* ::Input:: *)
(*addHalfULPInterval[nonnegative[aa]]*)


(* ::Input:: *)
(*\!\(\**)
(*TagBox[*)
(*RowBox[{"addHalfULPInterval", "[", *)
(*SuperscriptBox[*)
(*RowBox[{"nonnegative", "[", *)
(*RowBox[{*)
(*RowBox[{"Interval", "[", *)
(*RowBox[{"{", *)
(*RowBox[{"0", ",", "9"}], "}"}], "]"}], "+", *)
(*RowBox[{"Interval", "[", *)
(*RowBox[{"{", *)
(*RowBox[{*)
(*RowBox[{"-", *)
(*RowBox[{"Min", "[", *)
(*RowBox[{"0", ",", *)
(*FractionBox[*)
(*SuperscriptBox["2", *)
(*RowBox[{"Floor", "[", *)
(*FractionBox[*)
(*RowBox[{"Log", "[", "9", "]"}], *)
(*RowBox[{"Log", "[", "2", "]"}]], "]"}]], "9007199254740992"]}], "]"}]}], ",", *)
(*FractionBox[*)
(*SuperscriptBox["2", *)
(*RowBox[{"Floor", "[", *)
(*FractionBox[*)
(*RowBox[{"Log", "[", "9", "]"}], *)
(*RowBox[{"Log", "[", "2", "]"}]], "]"}]], "9007199254740992"]}], "}"}], "]"}]}], "]"}], "2"], "]"}],*)
(*HoldForm]\)//Trace*)


(* ::Input:: *)
(*Interval[{0,9}]+Interval[{-Min[0,2^Floor[Log[9]/Log[2]]/9007199254740992],2^Floor[Log[9]/Log[2]]/9007199254740992}]*)


(* ::Input:: *)
(*??\!\(\**)
(*TagBox["addHalfULPInterval",*)
(*HoldForm]\)*)


(* ::Text:: *)
(*The error analysis from [SZ05], section 2.1 (section 3.1 in the preprint, which is more explicit when it comes to the analysis.)*)


(* ::Input::Initialization:: *)
Begin["`SZ05`"]


(* ::Input::Initialization:: *)
bits[x_Interval]:=N[Log2[Max[Abs[x]]]];


(* ::Input::Initialization:: *)
a2=1/2


(* ::Input::Initialization:: *)
a3=CorrectlyRound[1/3!]


(* ::Input::Initialization:: *)
a4=CorrectlyRound[1/4!]


(* ::Input::Initialization:: *)
a5=CorrectlyRound[1/5!]


(* ::Input::Initialization:: *)
h=Interval[{CorrectlyRound[-(2^-10.977``100)],CorrectlyRound[2^-10.977``100]}]


(* ::Input::Initialization:: *)
k=IEEEEvaluateWithAbsoluteError[h^2]


(* ::Input::Initialization:: *)
\[Delta]k=Abs[IEEEEvaluateWithAbsoluteError[hh^2]-hh^2/.hh->h]


(* ::Input::Initialization:: *)
bits[k]


(* ::Input::Initialization:: *)
bits[\[Delta]k]


(* ::Input::Initialization:: *)
k'=IEEEEvaluateWithAbsoluteError[h^3]


(* ::Input::Initialization:: *)
\[Delta]k'=Abs[Collect[IEEEEvaluateWithAbsoluteError[hh^3]-hh^3,hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[k']


(* ::Input::Initialization:: *)
bits[\[Delta]k']


(* ::Text:: *)
(*TODO(phl): Propagate positiveness of the lower bound.*)


(* ::Input::Initialization:: *)
S\[FivePointedStar]1=IEEEEvaluateWithAbsoluteError[a5 k]


(* ::Input::Initialization:: *)
\[Delta]S\[FivePointedStar]1=Abs[Collect[IEEEEvaluateWithAbsoluteError[a5 hh^2]-a5 hh^2,hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]1]


(* ::Input::Initialization:: *)
bits[\[Delta]S\[FivePointedStar]1]


(* ::Input::Initialization:: *)
S\[FivePointedStar]2=IEEEEvaluateWithAbsoluteError[S\[FivePointedStar]1-a3]


(* ::Input::Initialization:: *)
\[Delta]S\[FivePointedStar]2=Abs[Collect[IEEEEvaluateWithAbsoluteError[a5 hh^2-a3]-(a5 hh^2-a3),hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]2]


(* ::Input::Initialization:: *)
bits[\[Delta]S\[FivePointedStar]2]


(* ::Input::Initialization:: *)
S\[FivePointedStar]3=IEEEEvaluateWithAbsoluteError[k' S\[FivePointedStar]2]


(* ::Input::Initialization:: *)
\[Delta]S\[FivePointedStar]3=Abs[Collect[IEEEEvaluateWithAbsoluteError[hh^3(a5 hh^2-a3)]-hh^3(a5 hh^2-a3),hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]3]


(* ::Input::Initialization:: *)
bits[\[Delta]S\[FivePointedStar]3]


(* ::Input::Initialization:: *)
C\[FivePointedStar]1=IEEEEvaluateWithAbsoluteError[a4 k]


(* ::Input::Initialization:: *)
\[Delta]C\[FivePointedStar]1=Abs[Collect[IEEEEvaluateWithAbsoluteError[a4 hh^2]-a4 hh^2,hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]1]


(* ::Input::Initialization:: *)
bits[\[Delta]C\[FivePointedStar]1]


(* ::Input::Initialization:: *)
C\[FivePointedStar]2=IEEEEvaluateWithAbsoluteError[C\[FivePointedStar]1-a2]


(* ::Text:: *)
(*TODO(phl): Because our intervals "bleed" a bit below 0 for values that should remain positive, we don't have the expected bound here and we lose at least a bit.  Error evaluation differ after this point.*)


(* ::Input::Initialization:: *)
Min[C\[FivePointedStar]2]>=-1/2


(* ::Input::Initialization:: *)
\[Delta]C\[FivePointedStar]2=Abs[Collect[IEEEEvaluateWithAbsoluteError[a4 hh^2-a2]-(a4 hh^2-a2),hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]2]


(* ::Input::Initialization:: *)
bits[\[Delta]C\[FivePointedStar]2]


(* ::Input::Initialization:: *)
C\[FivePointedStar]3=IEEEEvaluateWithAbsoluteError[k C\[FivePointedStar]2]


(* ::Input::Initialization:: *)
\[Delta]C\[FivePointedStar]3=Abs[Collect[IEEEEvaluateWithAbsoluteError[hh^2(a4 hh^2-a2)]-hh^2(a4 hh^2-a2),hh]/.hh->h]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]3]


(* ::Input::Initialization:: *)
bits[\[Delta]C\[FivePointedStar]3]


(* ::Input::Initialization:: *)
End[]
