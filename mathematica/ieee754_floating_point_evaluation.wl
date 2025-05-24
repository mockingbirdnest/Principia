(* ::Package:: *)

BeginPackage["IEEE754FloatingPointEvaluation`"]


Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]]


IEEEEvaluate;
IEEEEvaluate::usage =
"IEEEEvaluate[\!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\)] " <>
"evaluates \!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\) " <>
"according to the rules of IEEE arithmetic.  It uses the format set by " <>
"SetFloatingPointFormat and the rounding mode set by SetRoundingMode.  " <>
"It does not support associativity, i.e., parentheses are required.";
IEEEEvaluate::argnum =
"IEEEEvaluate called with `1` arguments; 1 argument is expected.";
IEEEEvaluate::badarg =
"IEEEEvaluate must be called with Plus, Times, Power, or a number."
IEEEEvaluate::badass =
"IEEEEvaluate does not support associativity, expressions must be " <>
"parenthesized."


IEEEEvaluateInterval;
IEEEEvaluateInterval::badass =
"IEEEEvaluateInterval does not support associativity, expressions must be " <>
"parenthesized."


IEEEEvaluateAbsoluteInterval;
IEEEEvaluateAbsoluteInterval::badass =
"IEEEEvaluateAbsoluteInterval does not support associativity, expressions must be " <>
"parenthesized."


UseFMA;
UseFMA::usage =
"UseFMA is an option for IEEEEvaluate that specifies whether to use FMA " <>
"for expressions of the form \!\(\*
StyleBox[\"a\",\nFontSlant->\"Italic\"]\) * \!\(\*
StyleBox[\"b\",\nFontSlant->\"Italic\"]\) + \!\(\*
StyleBox[\"c\",\nFontSlant->\"Italic\"]\).";


Begin["`Private`"]


SetAttributes[IEEEEvaluate,HoldAll];
Options[IEEEEvaluate]={UseFMA->True};
IEEEEvaluate[x:(_Plus|_Subtract|_Times|_Divide|_Power|_?NumberQ),OptionsPattern[]]:=
Block[
{Plus,Times,ev},
SetAttributes[ev,HoldAll];
ev[a_*b_+c_]:=
If[
OptionValue[UseFMA],
CorrectlyRound[IEEEEvaluate[a]IEEEEvaluate[b]+IEEEEvaluate[c]],
CorrectlyRound[CorrectlyRound[IEEEEvaluate[a]IEEEEvaluate[b]]+IEEEEvaluate[c]]];
ev[a_+b_]:=CorrectlyRound[IEEEEvaluate[a]+IEEEEvaluate[b]];
ev[a_+b__]:=(Message[IEEEEvaluate::badass]; $Failed);
ev[a_-b_]:=CorrectlyRound[IEEEEvaluate[a]-IEEEEvaluate[b]];
ev[a_*b_]:=CorrectlyRound[IEEEEvaluate[a]*IEEEEvaluate[b]];
ev[a_*b__]:=(Message[IEEEEvaluate::badass]; $Failed);
ev[a_/b_]:=CorrectlyRound[IEEEEvaluate[a]/IEEEEvaluate[b]];
ev[a_^2]:=CorrectlyRound[IEEEEvaluate[a ]IEEEEvaluate[a]];
ev[a_^3]:=CorrectlyRound[IEEEEvaluate[a^2]IEEEEvaluate[a]];
ev[a_^4]:=CorrectlyRound[IEEEEvaluate[a^2]IEEEEvaluate[a^2]];
ev[a_?NumberQ]:=CorrectlyRound[a];
ev[x]];
IEEEEvaluate[_]:=
(Message[IEEEEvaluate::badarg]; $Failed);
IEEEEvaluate[_, args__]:=
(Message[IEEEEvaluate::argnum, Length[{args}] + 1]; $Failed);


halfULPBelow1:=1-FromRepresentation[Representation[1]-1/2];


intervalMax[x_Interval]:=Max[x];
intervalMin[x_Interval]:=Min[x];


(*Would want to write Max[Log[...]] below but that doesn't work somehow because evai is HoldAll.*)
halfULP[x_]:=Block[{exponent=Log2[intervalMax[Abs[x]]]},halfULPBelow1 2^Ceiling[exponent]];


Options[addHalfULPInterval]={Positive->False};
addHalfULPInterval[x_,OptionsPattern[]]:=
Block[{h=halfULP[x]},If[OptionValue[Positive],x+Interval[{0,h}],x+Interval[{-h,h}]]];


SetAttributes[IEEEEvaluateAbsoluteInterval,HoldAll];
Options[IEEEEvaluateAbsoluteInterval]={UseFMA->True};
IEEEEvaluateAbsoluteInterval[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evai,usefma=OptionValue[UseFMA]},
SetAttributes[evai,HoldAll];
evai[a_*b_+c_]:=If[
usefma,
addHalfULPInterval[evai[a]evai[b]+evai[c]],
addHalfULPInterval[addHalfULPInterval[evai[a]evai[b]]+evai[c]]];
evai[a_+b_]:=addHalfULPInterval[evai[a]+evai[b]];
evai[a_+b__]:=(Message[IEEEEvaluateAbsoluteInterval::badass]; $Failed);
evai[a_*b_]:=addHalfULPInterval[evai[a]evai[b]];
evai[a_*b__]:=(Message[IEEEEvaluateAbsoluteInterval::badass]; $Failed);
evai[a_/b_]:=addHalfULPInterval[evai[a]/evai[b]];
(*Negation is exact.*)
evai[-a_]:=-evai[a];
(*Squaring an interval is not the same as multiplying two identical intervals.
Also, if the lower bound of the square is 0, it is exact.*) 
evai[a_^2]:=addHalfULPInterval[evai[a]^2,Positive->True];
evai[a_^3]:=Block[{t},addHalfULPInterval[addHalfULPInterval[evai[t]^2]evai[t]]/.t->a];
evai[a_^4]:=Block[
{t},
addHalfULPInterval[addHalfULPInterval[evai[t]^2,Positive->True]^2,Positive->True]/.t->a];
evai[a_?NumberQ]:=CorrectlyRound[a];
evai[a_]:=ReleaseHold[a];
evai[x]];


End[]


EndPackage[]
