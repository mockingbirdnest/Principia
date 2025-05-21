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


\[Delta]Interval:=Block[{abs\[Delta]=FromRepresentation[Representation[1]-1]-1},Interval[{-abs\[Delta],abs\[Delta]}]];


SetAttributes[IEEEEvaluateInterval,HoldAll];
Options[IEEEEvaluateInterval]={UseFMA->True};
IEEEEvaluateInterval[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evi,usefma=OptionValue[UseFMA]},
SetAttributes[evi,HoldAll];
evi[a_*b_+c_]:=If[usefma,(evi[a]evi[b]+evi[c])(1+\[Delta]Interval),((evi[a]evi[b])(1+\[Delta]Interval)+evi[c])(1+\[Delta]Interval)];
evi[a_+b_]:=(evi[a]+evi[b])(1+\[Delta]Interval);
evi[a_+b__]:=(Message[IEEEEvaluateInterval::badass]; $Failed);
evi[a_*b_]:=(evi[a]evi[b])(1+\[Delta]Interval);
evi[a_*b__]:=(Message[IEEEEvaluateInterval::badass]; $Failed);
evi[a_/b_]:=(evi[a]/evi[b])(1+\[Delta]Interval);
(*Squaring an interval is not the same as multiplying two identical intervals.*) 
evi[a_^2]:=evi[a]^2(1+\[Delta]Interval);
evi[a_^3]:=evi[a]^3(1+\[Delta]Interval)(1+\[Delta]Interval);
evi[a_^4]:=evi[a]^4(1+\[Delta]Interval)(1+\[Delta]Interval);
evi[a_?NumberQ]:=Block[{cra=CorrectlyRound[a]},Interval[{cra,cra}]];
evi[a_]:=ReleaseHold[a];
evi[x]];


intervalMax[x_Interval]:=Max[x];


halfULP[x_]:=Block[{
(*Would want to write Max[Log[...]] below but that doesn't work somehow because evai is HoldAll.*)
exponent=Log2[intervalMax[Abs[x]]],
halfULP1=1-FromRepresentation[Representation[1]-1/2]},
halfULP1 2^Ceiling[exponent]];


addHalfULPInterval[x_]:=Block[{h=halfULP[x]},x+Interval[{-h,h}]];


SetAttributes[IEEEEvaluateAbsoluteInterval,HoldAll];
Options[IEEEEvaluateAbsoluteInterval]={UseFMA->True};
IEEEEvaluateAbsoluteInterval[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evai,usefma=OptionValue[UseFMA]},
SetAttributes[evai,HoldAll];
evai[a_*b_+c_]:=If[usefma,(evai[a]evai[b]+evai[c])(1+\[Delta]Interval),((evai[a]evai[b])(1+\[Delta]Interval)+evai[c])(1+\[Delta]Interval)];
evai[a_+b_]:=addHalfULPInterval[evai[a]+evai[b]];
evai[a_+b__]:=(Message[IEEEEvaluateInterval::badass]; $Failed);
evai[a_*b_]:=addHalfULPInterval[evai[a]evai[b]];
evai[a_*b__]:=(Message[IEEEEvaluateInterval::badass]; $Failed);
evai[a_/b_]:=addHalfULPInterval[evai[a]/evai[b]];
(*Squaring an interval is not the same as multiplying two identical intervals.*) 
evai[a_^2]:=addHalfULPInterval[evai[a]^2];
evai[a_^3]:=Nest[addHalfULPInterval,evai[a]^3,2];
evai[a_^4]:=Nest[addHalfULPInterval,evai[a]^4,2];
evai[a_?NumberQ]:=Block[{cra=CorrectlyRound[a]},Interval[{cra,cra}]];
evai[a_]:=ReleaseHold[a];
evai[x]];


End[]


EndPackage[]
