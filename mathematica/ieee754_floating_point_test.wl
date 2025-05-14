(* ::Package:: *)

(* ::Input::Initialization:: *)
Get[FileNameJoin[NotebookDirectory[],"ieee754_floating_point.wl"]]


(* ::Input::Initialization:: *)
On[Assert]


\[LeftAngleBracket]x_\[RightAngleBracket]:=CorrectlyRound[x];


SetFloatingPointFormat[binary32];
SetRoundingMode[NearestTiesToEven];
inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];


UlpPlot[{
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},
PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary32, rounding to nearest, ties to even",
ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]


(* ::Input::Initialization:: *)
Assert[Bits[Sqrt[2]]=="0|01111111|01101010000010011110011;0011001111\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]]=="0|01111111|01101010000010011110011;"]
Assert[Bits[-Sqrt[2]]=="1|01111111|01101010000010011110011;0011001111\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]=="1|01111111|01101010000010011110011;"]


SetFloatingPointFormat[binary32];
SetRoundingMode[TowardPositiveInfinity]


UlpPlot[{UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},
PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary32, rounding toward +\[Infinity]",
ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]


(* ::Input::Initialization:: *)
Assert[Bits[Sqrt[2]]=="0|01111111|01101010000010011110011;0011001111\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]]=="0|01111111|01101010000010011110100;"]
Assert[Bits[-Sqrt[2]]=="1|01111111|01101010000010011110011;0011001111\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]=="1|01111111|01101010000010011110100;"]


SetFloatingPointFormat[binary64];
SetRoundingMode[NearestTiesToEven];
inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];


UlpPlot[{UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,
UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},
PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary64, rounding to nearest, ties to even",
ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]


(* ::Input::Initialization:: *)
Assert[Bits[Sqrt[2]]=="0|01111111111|0110101000001001111001100110011111110011101111001100;1001000010\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]]=="0|01111111111|0110101000001001111001100110011111110011101111001101;"]
Assert[Bits[-Sqrt[2]]=="1|01111111111|0110101000001001111001100110011111110011101111001100;1001000010\[Ellipsis]"]
Assert[Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]=="1|01111111111|0110101000001001111001100110011111110011101111001101;"]


Table[
SetFloatingPointFormat[format];
{#,{Bits[#],HexLiteral[#]},{Bits[CorrectlyRound[#]],HexLiteral[CorrectlyRound[#]]}}&/@{\[Pi],1/3,Sqrt[2],1+2^-format[[1]],1+2^(-format[[1]]-5)}//TableForm,
{format,{binary16,binary32,binary64,x87extended}}]//TableForm


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary64];
SetRoundingMode[TowardPositiveInfinity];
Assert[CorrectlyRound[1/3]==3002399751580331/9007199254740992]
SetRoundingMode[TowardNegativeInfinity];
Assert[CorrectlyRound[1/3]==6004799503160661/18014398509481984]
Assert[CorrectlyRound[1/3,RoundingMode->Toward0]==6004799503160661/18014398509481984]
Assert[CorrectlyRound[1/3,RoundingMode->TowardPositiveInfinity]==3002399751580331/9007199254740992]
Assert[CorrectlyRound[1/3,RoundingMode->TowardNegativeInfinity]==6004799503160661/18014398509481984]
Assert[CorrectlyRound[1/3,RoundingMode->NearestTiesToEven]==6004799503160661/18014398509481984]
