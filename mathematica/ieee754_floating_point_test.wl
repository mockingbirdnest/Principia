(* ::Package:: *)

(* ::Input::Initialization:: *)
Get[FileNameJoin[NotebookDirectory[],"ieee754_floating_point.wl"]]


(* ::Input::Initialization:: *)
On[Assert]


(* ::Input::Initialization:: *)
\[LeftAngleBracket]x_\[RightAngleBracket]:=CorrectlyRound[x];


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary32];
SetRoundingMode[NearestTiesToEven];
inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];


(* ::Input::Initialization:: *)
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


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary32];
SetRoundingMode[TowardPositiveInfinity]


(* ::Input::Initialization:: *)
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


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary64];
SetRoundingMode[NearestTiesToEven];
inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];


(* ::Input::Initialization:: *)
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


(* ::Input::Initialization:: *)
Assert[Bits[Truncate[0,Sqrt[2]]]=="0|01111111111|0110101000001001111001100110011111110011101111001100;"]
Assert[Bits[Truncate[1,Sqrt[2]]]=="0|01111111111|0110101000001001111001100110011111110011101111001100;"]
Assert[Bits[Truncate[4,Sqrt[2]]]=="0|01111111111|0110101000001001111001100110011111110011101111000000;"]
Assert[Bits[Truncate[10,-Sqrt[2]]]=="1|01111111111|0110101000001001111001100110011111110011100000000000;"]


(* ::Input::Initialization:: *)
checkFormat[format_]:=(
SetFloatingPointFormat[format];SetRoundingMode[NearestTiesToEven];
{#,{Bits[#],HexLiteral[#],HexLiteral[#,Quotes->4]},{Bits[CorrectlyRound[#]],HexLiteral[CorrectlyRound[#]],HexLiteral[CorrectlyRound[#],Quotes->5]}}&/@
{\[Pi],1/3,Sqrt[2],1+2^-format[[1]],1+2^(-format[[1]]-5)})


(* ::Input::Initialization:: *)
Assert[checkFormat[binary16]==
{{Pi,{"0|10000|1001001000;0111111011\[Ellipsis]","0x1.92'1'FB5\[Ellipsis]p1_11_sigbits","0x1.921F'B5\[Ellipsis]p1_11_sigbits"},{"0|10000|1001001000;","0x1.92'0p1_11_sigbits","0x1.920p1_11_sigbits"}},{1/3,{"0|01101|0101010101;0101010101\[Ellipsis]","0x1.55'5'555\[Ellipsis]p-2_11_sigbits","0x1.5555'55\[Ellipsis]p-2_11_sigbits"},{"0|01101|0101010101;","0x1.55'4p-2_11_sigbits","0x1.554p-2_11_sigbits"}},{Sqrt[2],{"0|01111|0110101000;0010011110\[Ellipsis]","0x1.6A'0'9E6\[Ellipsis]p0_11_sigbits","0x1.6A09'E6\[Ellipsis]p0_11_sigbits"},{"0|01111|0110101000;","0x1.6A'0p0_11_sigbits","0x1.6A0p0_11_sigbits"}},{2049/2048,{"0|01111|0000000000;1000000000\[Ellipsis]","0x1.00'2'000\[Ellipsis]p0_11_sigbits","0x1.0020'00\[Ellipsis]p0_11_sigbits"},{"0|01111|0000000000;","0x1.00'0p0_11_sigbits","0x1.000p0_11_sigbits"}},{65537/65536,{"0|01111|0000000000;0000010000\[Ellipsis]","0x1.00'0'100\[Ellipsis]p0_11_sigbits","0x1.0001'00\[Ellipsis]p0_11_sigbits"},{"0|01111|0000000000;","0x1.00'0p0_11_sigbits","0x1.000p0_11_sigbits"}}}]


(* ::Input::Initialization:: *)
Assert[checkFormat[binary32]==
{{Pi,{"0|10000000|10010010000111111011010;1010001000\[Ellipsis]","0x1.921FB'5'444\[Ellipsis]p1f","0x1.921F'B544'4\[Ellipsis]p1f"},{"0|10000000|10010010000111111011011;","0x1.921FB'6p1f","0x1.921FB'6p1f"}},{1/3,{"0|01111101|01010101010101010101010;1010101010\[Ellipsis]","0x1.55555'5'555\[Ellipsis]p-2f","0x1.5555'5555'5\[Ellipsis]p-2f"},{"0|01111101|01010101010101010101011;","0x1.55555'6p-2f","0x1.55555'6p-2f"}},{Sqrt[2],{"0|01111111|01101010000010011110011;0011001111\[Ellipsis]","0x1.6A09E'6'67F\[Ellipsis]p0f","0x1.6A09'E667'F\[Ellipsis]p0f"},{"0|01111111|01101010000010011110011;","0x1.6A09E'6p0f","0x1.6A09E'6p0f"}},{16777217/16777216,{"0|01111111|00000000000000000000000;1000000000\[Ellipsis]","0x1.00000'1'000\[Ellipsis]p0f","0x1.0000'0100'0\[Ellipsis]p0f"},{"0|01111111|00000000000000000000000;","0x1.00000'0p0f","0x1.00000'0p0f"}},{536870913/536870912,{"0|01111111|00000000000000000000000;0000010000\[Ellipsis]","0x1.00000'0'080\[Ellipsis]p0f","0x1.0000'0008'0\[Ellipsis]p0f"},{"0|01111111|00000000000000000000000;","0x1.00000'0p0f","0x1.00000'0p0f"}}}]


(* ::Input::Initialization:: *)
Assert[checkFormat[binary64]==
{{Pi,{"0|10000000000|1001001000011111101101010100010001000010110100011000;0100011010\[Ellipsis]","0x1.921FB54442D18'469\[Ellipsis]p1","0x1.921F'B544'42D1'8469\[Ellipsis]p1"},{"0|10000000000|1001001000011111101101010100010001000010110100011000;","0x1.921FB54442D18p1","0x1.921FB'54442'D18p1"}},{1/3,{"0|01111111101|0101010101010101010101010101010101010101010101010101;0101010101\[Ellipsis]","0x1.5555555555555'555\[Ellipsis]p-2","0x1.5555'5555'5555'5555\[Ellipsis]p-2"},{"0|01111111101|0101010101010101010101010101010101010101010101010101;","0x1.5555555555555p-2","0x1.55555'55555'555p-2"}},{Sqrt[2],{"0|01111111111|0110101000001001111001100110011111110011101111001100;1001000010\[Ellipsis]","0x1.6A09E667F3BCC'908\[Ellipsis]p0","0x1.6A09'E667'F3BC'C908\[Ellipsis]p0"},{"0|01111111111|0110101000001001111001100110011111110011101111001101;","0x1.6A09E667F3BCDp0","0x1.6A09E'667F3'BCDp0"}},{9007199254740993/9007199254740992,{"0|01111111111|0000000000000000000000000000000000000000000000000000;1000000000\[Ellipsis]","0x1.0000000000000'800\[Ellipsis]p0","0x1.0000'0000'0000'0800\[Ellipsis]p0"},{"0|01111111111|0000000000000000000000000000000000000000000000000000;","0x1.0000000000000p0","0x1.00000'00000'000p0"}},{288230376151711745/288230376151711744,{"0|01111111111|0000000000000000000000000000000000000000000000000000;0000010000\[Ellipsis]","0x1.0000000000000'040\[Ellipsis]p0","0x1.0000'0000'0000'0040\[Ellipsis]p0"},{"0|01111111111|0000000000000000000000000000000000000000000000000000;","0x1.0000000000000p0","0x1.00000'00000'000p0"}}}]


(* ::Input::Initialization:: *)
Assert[checkFormat[x87extended]==
{{Pi,{"0|100000000000000|100100100001111110110101010001000100001011010001100001000110100;1100010011\[Ellipsis]","0x1.921FB54442D1846'9'898\[Ellipsis]p1l","0x1.921F'B544'42D1'8469'898\[Ellipsis]p1l"},{"0|100000000000000|100100100001111110110101010001000100001011010001100001000110101;","0x1.921FB54442D1846'Ap1l","0x1.921FB'54442'D1846'Ap1l"}},{1/3,{"0|011111111111101|010101010101010101010101010101010101010101010101010101010101010;1010101010\[Ellipsis]","0x1.555555555555555'5'555\[Ellipsis]p-2l","0x1.5555'5555'5555'5555'555\[Ellipsis]p-2l"},{"0|011111111111101|010101010101010101010101010101010101010101010101010101010101011;","0x1.555555555555555'6p-2l","0x1.55555'55555'55555'6p-2l"}},{Sqrt[2],{"0|011111111111111|011010100000100111100110011001111111001110111100110010010000100;0101100101\[Ellipsis]","0x1.6A09E667F3BCC90'8'B2F\[Ellipsis]p0l","0x1.6A09'E667'F3BC'C908'B2F\[Ellipsis]p0l"},{"0|011111111111111|011010100000100111100110011001111111001110111100110010010000100;","0x1.6A09E667F3BCC90'8p0l","0x1.6A09E'667F3'BCC90'8p0l"}},{18446744073709551617/18446744073709551616,{"0|011111111111111|000000000000000000000000000000000000000000000000000000000000000;1000000000\[Ellipsis]","0x1.000000000000000'1'000\[Ellipsis]p0l","0x1.0000'0000'0000'0001'000\[Ellipsis]p0l"},{"0|011111111111111|000000000000000000000000000000000000000000000000000000000000000;","0x1.000000000000000'0p0l","0x1.00000'00000'00000'0p0l"}},{590295810358705651713/590295810358705651712,{"0|011111111111111|000000000000000000000000000000000000000000000000000000000000000;0000010000\[Ellipsis]","0x1.000000000000000'0'080\[Ellipsis]p0l","0x1.0000'0000'0000'0000'080\[Ellipsis]p0l"},{"0|011111111111111|000000000000000000000000000000000000000000000000000000000000000;","0x1.000000000000000'0p0l","0x1.00000'00000'00000'0p0l"}}}]


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
