(* ::Package:: *)

BeginPackage["IEEE754FloatingPoint`"]


(*Useful, but doesn't really belong here...*)
ConditionNumber[f_,params_]:=D[f,{params}]params/f;


SetFloatingPointFormat;
binary256={237,19};
binary128={113,15};
(*|x87extended| should have the right arithmetic, but |Bits| will not show the demented explicit bit.*)
x87extended={64,15};
binary64={53,11};
binary32={24,8};
binary16={11,5};


SetRoundingMode;
NearestTiesToEven;
Toward0;
TowardPositiveInfinity;
TowardNegativeInfinity;


Representation;
FromRepresentation;


UlpDistance;
UlpPlot;


Bits;


HexLiteral;


CorrectlyRound;
\[LeftAngleBracket]x_\[RightAngleBracket]:=CorrectlyRound[x];


Begin["`Private`"]


SetFloatingPointFormat[{sBits_,eBits_}]:=
(significandBits=sBits;
 exponentBits=eBits;
 bias=2^(exponentBits-1)-1;)


SetRoundingMode[mode_]:=(correctlyRoundRepresentation=mode;);
NearestTiesToEven=Round;
Toward0=IntegerPart;
TowardPositiveInfinity=Ceiling;
TowardNegativeInfinity=Floor;


mantissaExponent12[x_]:={#[[1]]*2,#[[2]]-1}&[MantissaExponent[x,2]];
Representation[x_]:=Block[
{sign,magnitude,\[Mu],e},
sign=Sign[x];
magnitude=Abs[x];
If[magnitude==\[Infinity],
sign(2^exponentBits-1)2^(significandBits-1),
If[x==0,0,
{\[Mu],e}=mantissaExponent12[magnitude];
If[e<=-bias,
2^(significandBits-1) 2^(e+bias-1) \[Mu],
sign(2^(significandBits-1) (\[Mu]-1)+2^(significandBits-1) (e+bias))
]]]];
FromRepresentation[n_]:=Block[
{sign,magnitude,\[Mu],e},
sign=Sign[n];
magnitude=Abs[n];
\[Mu]=Mod[magnitude,2^(significandBits-1)];
e=IntegerPart[magnitude/2^(significandBits-1)];
If[e==0,
\[Mu]/2^(significandBits-1) 2^(1-bias),
sign(1+\[Mu]/2^(significandBits-1))2^(e-bias)]];
CorrectlyRound[x_]:=If[x==\[Infinity]||x==-\[Infinity],x,If[Abs[#]>=2^(bias+1),Sign[x]\[Infinity],#]&@FromRepresentation[correctlyRoundRepresentation[Representation[x]]]];
UlpDistance[x_,y_]:=Abs[Representation[x]-Representation[y]]


Bits[n_]:=If[n>=0,"0","1"]<>
"|"<>IntegerString[IntegerPart[Representation[Abs[n]]/2^(significandBits-1)],2,exponentBits]<>
"|"<>IntegerString[Mod[IntegerPart[Representation[Abs[n]]],2^(significandBits-1)],2,significandBits-1]<>
";"<>If[FractionalPart[Representation[Abs[n]]]==0,"",ToString/@RealDigits[N[FractionalPart[Representation[Abs[n]]],5],2,10,-1][[1]]<>"\[Ellipsis]"];


fullHexDigits:=Floor[(significandBits-1)/4]


leastFullHexDigitValue:=2^(significandBits-1)/16^fullHexDigits


leastHexDigitValue:=If[leastFullHexDigitValue>1,leastFullHexDigitValue/16,leastFullHexDigitValue]


HexLiteral[n_]:=If[n<0,"-",""]<>
"0x1."<>ToUpperCase[
IntegerString[
Mod[IntegerPart[Representation[Abs[n]]/leastFullHexDigitValue],16^fullHexDigits],16,fullHexDigits]<>
If[
leastHexDigitValue<1,
"'"<>IntegerString[Mod[IntegerPart[Representation[Abs[n]]/leastHexDigitValue],16],16,1],
""]<>
If[FractionalPart[Representation[Abs[n]]]==0,
"",
"'"<>ToString/@RealDigits[
N[FractionalPart[Representation[Abs[n]]/leastHexDigitValue],5],
16,3,-1][[1]]<>"\[Ellipsis]"]]<>
"p"<>ToString[IntegerPart[Representation[Abs[n]]/2^(significandBits-1)]-bias]<>
Switch[
{significandBits,exponentBits},
binary32,"f",
binary64,"",
x87extended,"l",
_,"_"<>ToString[significandBits]<>"_sigbits"]


smol=-12;
Log2OrSmol[x_]:=If[x==0,smol,If[Log2[x]<smol+3,smol+1,Log2[x]]];
UlpPlot[errors_,options___]:=Histogram[
Map[Log2OrSmol,N[errors,10],{2}],
{Join[{smol,smol+1/8},Table[k,{k,smol+3,significandBits+1,1/8}]]},
"Intensity",
Ticks->{Table[{k,If[k==smol,0,If[k==smol+1||k==smol+2,"",If[2^k>=1&&2^k<100,2^k,("2")^ToString[k]]]]},{k,-100,100}],{}},
ImageSize->1400,AspectRatio->1/4,PlotRange->{{smol,significandBits},Full},AxesLabel->{"ULPs",""},options]


End[]


EndPackage[]


(* ::Section:: *)
(*Examples*)


(* ::Code:: *)
(*SetFloatingPointFormat[binary32];*)
(*SetRoundingMode[NearestTiesToEven];*)
(*inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];*)


(* ::Code:: *)
(*UlpPlot[{*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},*)
(*PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary32, rounding to nearest, ties to even",*)
(*ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]*)


(* ::Code:: *)
(*{Bits[Sqrt[2]],Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]],Bits[-Sqrt[2]],Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]}//Column*)


(* ::Code:: *)
(*SetFloatingPointFormat[binary32];*)
(*SetRoundingMode[TowardPositiveInfinity]*)


(* ::Code:: *)
(*UlpPlot[{UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},*)
(*PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary32, rounding toward +\[Infinity]",*)
(*ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]*)


(* ::Code:: *)
(*{Bits[Sqrt[2]],Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]],Bits[-Sqrt[2]],Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]}//Column*)


(* ::Code:: *)
(*SetFloatingPointFormat[binary64];*)
(*SetRoundingMode[NearestTiesToEven];*)
(*inputs=Reverse/@Sort/@Transpose[{CorrectlyRound/@Range[0,1,1/1000],CorrectlyRound/@Range[1/1*^7,1+1/1*^7,1/1000]}];*)


(* ::Code:: *)
(*UlpPlot[{UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]-#[[2]]\[RightAngleBracket]\[LeftAngleBracket]#[[1]]+#[[2]]\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[\[LeftAngleBracket]\[LeftAngleBracket]#[[1]]^2\[RightAngleBracket]-\[LeftAngleBracket]#[[2]]^2\[RightAngleBracket]\[RightAngleBracket]]\[RightAngleBracket]]&/@inputs,*)
(*UlpDistance[Sqrt[#[[1]]^2-#[[2]]^2],\[LeftAngleBracket]Sqrt[#[[1]]^2-#[[2]]^2]\[RightAngleBracket]]&/@inputs},*)
(*PlotLabel->"Forward error in the computation of \!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\) in binary64, rounding to nearest, ties to even",*)
(*ChartLegends->{"\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]x - y\[RightAngleBracket] \[LeftAngleBracket]x + y\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\[LeftAngleBracket]\[LeftAngleBracket]\*SuperscriptBox[\(x\), \(2\)]\[RightAngleBracket] - \[LeftAngleBracket]\*SuperscriptBox[\(y\), \(2\)]\[RightAngleBracket]\[RightAngleBracket]\)]\)\[RightAngleBracket]","\[LeftAngleBracket]\!\(\*SqrtBox[\(\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)]\)]\)\[RightAngleBracket]"}]*)


(* ::Code:: *)
(*{Bits[Sqrt[2]],Bits[\[LeftAngleBracket]Sqrt[2]\[RightAngleBracket]],Bits[-Sqrt[2]],Bits[\[LeftAngleBracket]-Sqrt[2]\[RightAngleBracket]]}//Column*)


(* ::Code:: *)
(*Table[*)
(*SetFloatingPointFormat[format];*)
(*{#,{Bits[#],HexLiteral[#]},{Bits[CorrectlyRound[#]],HexLiteral[CorrectlyRound[#]]}}&/@{\[Pi],1/3,Sqrt[2],1+2^-format[[1]],1+2^(-format[[1]]-5)}//TableForm,*)
(*{format,{binary16,binary32,binary64,x87extended}}]//TableForm*)
