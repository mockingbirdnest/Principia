(* ::Package:: *)

BeginPackage["IEEE754FloatingPoint`"]


(*Useful, but doesn't really belong here...*)
ConditionNumber;


SetFloatingPointFormat;
SetFloatingPointFormat::usage =
"SetFloatingPointFormat[\!\(\*StyleBox[\"format\", \"TI\"]\)] defines the " <>
"floating-point format to use; allowed values are " <>
"\!\(\*StyleBox[\"binary16\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"binary32\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"binary64\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"x87extended\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"binary128\", \"TI\"]\), and " <>
"\!\(\*StyleBox[\"binary256\", \"TI\"]\).";
SetFloatingPointFormat::argnum =
"SetFloatingPointFormat called with `1` arguments; 1 argument is expected.";

binary256={237,19};
binary128={113,15};
(*|x87extended| should have the right arithmetic, but |Bits| will not show the demented explicit bit.*)
x87extended={64,15};
binary64={53,11};
binary32={24,8};
binary16={11,5};


SetRoundingMode;
SetRoundingMode::usage =
"SetRoundingMode[\!\(\*StyleBox[\"mode\", \"TI\"]\)] defines the " <>
"floating-point rounding mode to use; allowed values are " <>
"\!\(\*StyleBox[\"NearestTiesToEven\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"Toward0\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"TowardPositiveInfinity\", \"TI\"]\), and " <>
"\!\(\*StyleBox[\"TowardNegativeInfinity\", \"TI\"]\).";
SetRoundingMode::argnum =
"SetRoundingMode called with `1` arguments; 1 argument is expected.";
NearestTiesToEven;
Toward0;
TowardPositiveInfinity;
TowardNegativeInfinity;


RoundingMode;
RoundingMode::usage =
"RoundingMode is an option for " <>
"\!\(\*StyleBox[\"CorrectlyRound\", \"TI\"]\) that specifies the rounding " <>
"mode to use; allowed values are " <>
"\!\(\*StyleBox[\"NearestTiesToEven\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"Toward0\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"TowardPositiveInfinity\", \"TI\"]\), and " <>
"\!\(\*StyleBox[\"TowardNegativeInfinity\", \"TI\"]\).";


Representation;
Representation::usage =
"Representation[\!\(\*StyleBox[\"x\", \"TI\"]\)] returns the " <>
"machine representation of \!\(\*StyleBox[\"x\", \"TI\"]\) in the format " <>
"specified by \!\(\*StyleBox[\"SetFloatingPointFormat\", \"TI\"]\).  The " <>
"integer part of the result corresponds to the IEEE " <>
"(\!\(\*StyleBox[\"sign\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"biased exponent\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"trailing significand field\", \"TI\"]\)) as an unsigned " <>
"integer, its fractional part to the bits after the significand."
Representation::argnum =
"Representation called with `1` arguments; 1 argument is expected.";


FromRepresentation;
FromRepresentation::usage =
"FromRepresentation[\!\(\*StyleBox[\"repr\", \"TI\"]\)] returns " <>
"the number whose representation is \!\(\*StyleBox[\"repr\", \"TI\"]\).  " <>
"\!\(\*StyleBox[\"repr\", \"TI\"]\)  must be an unsigned integer " <>
"corresponding to the IEEE (\!\(\*StyleBox[\"sign\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"biased exponent\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"trailing significand field\", \"TI\"]\)), possibly with a " <>
"fractional part, similar to the result of " <>
"\!\(\*StyleBox[\"FromRepresentation\", \"TI\"]\).";
FromRepresentation::argnum =
"FromRepresentation called with `1` arguments; 1 argument is expected.";


UlpDistance;
UlpDistance::usage =
"UlpDistance[\!\(\*StyleBox[\"x\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"y\", \"TI\"]\)] returns the distance between " <>
"\!\(\*StyleBox[\"x\", \"TI\"]\) and \!\(\*StyleBox[\"y\", \"TI\"]\), " <>
"expressed in (possibly fractional) ULPs.  It is a distance, and therefore " <>
"nonnegative."
UlpDistance::argnum =
"UlpDistance called with `1` arguments; 2 arguments are expected.";


UlpPlot;
UlpPlot::usage =
"UlpPlot[{\!\(\*StyleBox[\"errors\", \"TI\"]\)}, " <>
"\!\(\*StyleBox[\"options\", \"TI\"]\)...] produces a histogram of the " <>
"\!\(\*StyleBox[\"errors\", \"TI\"]\).  " <>
"The \!\(\*StyleBox[\"options\", \"TI\"]\) are passed to the Histogram " <>
"function."
UlpPlot::argnum =
"UlpPlot called with `1` arguments; at least 1 argument is expected.";


Bits;
Bits::usage =
"Bits[\!\(\*StyleBox[\"x\", \"TI\"]\)] or " <>
"Bits[\!\(\*StyleBox[\"x\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"extraBits\", \"TI\"]\)] returns a string representation " <>
"of the bits of \!\(\*StyleBox[\"x\", \"TI\"]\), separating out the " <>
"\!\(\*StyleBox[\"sign\", \"TI\"]\), " <>
"\!\(\*StyleBox[\"biased exponent\", \"TI\"]\), " <>
"and \!\(\*StyleBox[\"trailing significand field\", \"TI\"]\).  If " <>
"\!\(\*StyleBox[\"extraBits\", \"TI\"]\) is provided, it specifies the " <>
"number of bits to display after the significand; it defaults to 10."
Bits::argnum =
"Bits called with `1` arguments; 1 or 2 arguments are expected.";


HexLiteral;
HexLiteral::usage =
"HexLiteral[\!\(\*StyleBox[\"x\", \"TI\"]\)] returns a string " <>
"representation of \!\(\*StyleBox[\"x\", \"TI\"]\) suitable for use in the " <>
"C++ language.  The digits after the significand are separated by a single " <>
"tick."
HexLiteral::argnum =
"HexLiteral called with `1` arguments; 1 argument is expected.";
Attributes[HexLiteral]={Listable}


CorrectlyRound;
CorrectlyRound::usage =
"CorrectlyRound[\!\(\*StyleBox[\"x\", \"TI\"]\)] returns a signed rational " <>
"obtained by rounding \!\(\*StyleBox[\"x\", \"TI\"]\) according to the " <>
"rounding mode set by SetRoundingMode and the format set by " <>
"SetFloatingPointFormat, or by the option RoundingMode."
CorrectlyRound::argnum =
"CorrectlyRound called with `1` arguments; 1 argument is expected.";
Attributes[CorrectlyRound]={Listable};


Truncate;
Truncate::usage =
"Truncate[\!\(\*StyleBox[\"\[Kappa]\", \"TI\"]\),\!\(\*StyleBox[\"x\", \"TI\"]\)] " <>
"returns the value obtained by clearing the last \!\(\*StyleBox[\"\[Kappa]\", \"TI\"]\) " <>
"bits of the mantissa of \!\(\*StyleBox[\"x\", \"TI\"]\)."


Begin["`Private`"]


ConditionNumber[f_,params_]:=D[f,{params}]params/f;


SetFloatingPointFormat[{sBits_,eBits_}]:=
(significandBits=sBits;
 exponentBits=eBits;
 bias=2^(exponentBits-1)-1;)
SetFloatingPointFormat[args___] :=
(Message[SetFloatingPointFormat::argnum, Length[{args}]]; $Failed)



SetRoundingMode[mode_]:=(correctlyRoundRepresentation=mode;);
SetRoundingMode[args___] :=
(Message[SetRoundingMode::argnum, Length[{args}]]; $Failed)
NearestTiesToEven=Round;
Toward0=IntegerPart;
TowardPositiveInfinity=Ceiling;
TowardNegativeInfinity=Floor;


RoundingMode[] := correctlyRoundRepresentation;


mantissaExponent12[x_]:={#[[1]]*2,#[[2]]-1}&[MantissaExponent[x,2]];

Representation[x_]:=Block[
{sign,magnitude,\[Mu],e},
sign=If[x<0,2^(exponentBits+significandBits-1),0];
magnitude=Abs[x];
If[magnitude==\[Infinity],
sign+(2^exponentBits-1)2^(significandBits-1),
If[x==0,0,
{\[Mu],e}=mantissaExponent12[magnitude];
If[e<=-bias,
2^(significandBits-1) 2^(e+bias-1) \[Mu],
sign+(2^(significandBits-1) (\[Mu]-1)+2^(significandBits-1) (e+bias))
]]]];
Representation[args___] :=
(Message[Representation::argnum, Length[{args}]]; $Failed)


FromRepresentation[n_]:=Block[
{sign,magnitude,\[Mu],e},
sign=If[n>=2^(exponentBits+significandBits-1),-1,1];
magnitude=Mod[n,2^(exponentBits+significandBits-1)];
\[Mu]=Mod[magnitude,2^(significandBits-1)];
e=IntegerPart[magnitude/2^(significandBits-1)];
If[e==0,
\[Mu]/2^(significandBits-1) 2^(1-bias),
sign(1+\[Mu]/2^(significandBits-1))2^(e-bias)]];
FromRepresentation[args___] :=
(Message[FromRepresentation::argnum, Length[{args}]]; $Failed)


UlpDistance[x_,y_]:=Abs[Representation[x]-Representation[y]];
UlpDistance[args___] :=
(Message[UlpDistance::argnum, Length[{args}]]; $Failed)


smol=-12;
Log2OrSmol[x_]:=If[x==0,smol,If[Log2[x]<smol+3,smol+1,Log2[x]]];
UlpPlot[errors_,options___]:=Histogram[
Map[Log2OrSmol,N[errors,10],{2}],
{Join[{smol,smol+1/8},Table[k,{k,smol+3,significandBits+1,1/8}]]},
"Intensity",
Ticks->{Table[{k,If[k==smol,0,If[k==smol+1||k==smol+2,"",If[2^k>=1&&2^k<100,2^k,("2")^ToString[k]]]]},{k,-100,100}],{}},
ImageSize->1400,AspectRatio->1/4,PlotRange->{{smol,significandBits},Full},AxesLabel->{"ULPs",""},options]
UlpPlot[] :=
(Message[UlpPlot::argnum, 0]; $Failed)


Bits[x_, extraBits_: 10]:=If[x>=0,"0","1"]<>
"|"<>IntegerString[IntegerPart[Representation[Abs[x]]/2^(significandBits-1)],2,exponentBits]<>
"|"<>IntegerString[Mod[IntegerPart[Representation[Abs[x]]],2^(significandBits-1)],2,significandBits-1]<>
";"<>If[FractionalPart[Representation[Abs[x]]]==0,"",ToString/@RealDigits[N[FractionalPart[Representation[Abs[x]]],extraBits/2],2,extraBits,-1][[1]]<>"\[Ellipsis]"];
Bits[args___] :=
(Message[Bits::argnum, Length[{args}]]; $Failed)


fullHexDigits:=Floor[(significandBits-1)/4]
leastFullHexDigitValue:=2^(significandBits-1)/16^fullHexDigits
leastHexDigitValue:=If[leastFullHexDigitValue>1,leastFullHexDigitValue/16,leastFullHexDigitValue]
HexLiteral[x_]:=If[x==0,"0.0",If[x<0,"-",""]<>
"0x1."<>ToUpperCase[
IntegerString[
Mod[IntegerPart[Representation[Abs[x]]/leastFullHexDigitValue],16^fullHexDigits],16,fullHexDigits]<>
If[
leastHexDigitValue<1,
"'"<>IntegerString[Mod[IntegerPart[Representation[Abs[x]]/leastHexDigitValue],16],16,1],
""]<>
If[FractionalPart[Representation[Abs[x]]]==0,
"",
"'"<>(IntegerString[#,16]&)/@RealDigits[
N[FractionalPart[Representation[Abs[x]]/leastHexDigitValue],5],
16,3,-1][[1]]<>"\[Ellipsis]"]]<>
"p"<>ToString[IntegerPart[Representation[Abs[x]]/2^(significandBits-1)]-bias]<>
Switch[
{significandBits,exponentBits},
binary32,"f",
binary64,"",
x87extended,"l",
_,"_"<>ToString[significandBits]<>"_sigbits"]]
HexLiteral[args___] :=
(Message[HexLiteral::argnum, Length[{args}]]; $Failed)


CorrectlyRound[x_,OptionsPattern[]]:=With[
{rounding=If[OptionValue[RoundingMode]===Automatic,correctlyRoundRepresentation,OptionValue[RoundingMode]]},
If[x==\[Infinity]||x==-\[Infinity],x,
If[Abs[#]>=2^(bias+1),Sign[x]\[Infinity],#]&
@FromRepresentation[rounding[Representation[x]]]]];
Options[CorrectlyRound]={RoundingMode->Automatic};
CorrectlyRound[args___] :=
(Message[CorrectlyRound::argnum, Length[{args}]]; $Failed)


Truncate[\[Kappa]_,x_]:=FromRepresentation[Floor[Representation[x], 2^\[Kappa]]]


End[]


EndPackage[]
