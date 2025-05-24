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


IEEEEvaluateWithAbsoluteError;
IEEEEvaluateWithAbsoluteError::usage =
"IEEEEvaluateWithAbsoluteError[\!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\)] " <>
"evaluates \!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\), " <>
"which may include intervals or unbound variables, with proper " <>
"propagation of absolute error (1/2 ULP of the result on each operation)."
IEEEEvaluateWithAbsoluteError::argnum =
"IEEEEvaluate called with `1` arguments; 1 argument is expected.";
IEEEEvaluateWithAbsoluteError::badass =
"IEEEEvaluateWithAbsoluteError does not support associativity, expressions must be " <>
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


halfULPAbove1:=FromRepresentation[Representation[1]+1/2]-1;


(* ::Text:: *)
(*Min and Max when applied to an undefined value return that value.  These functions don't do that, they just stay unevaluated.*)


intervalMax[x_Interval]:=Max[x];
intervalMin[x_Interval]:=Min[x];


(* ::Text:: *)
(*Returns a half ULP above for the largest element (in absolute value) of its argument.  The returned value is positive, regardless of the sign of the argument.  The argument is an interval or an unbound variable.*)


(* Would want to write Max[Log[...]] below but that doesn't work somehow because evae is HoldAll. *)
halfULP[x_]:=Block[{exponent=Log2[intervalMax[Abs[x]]]},halfULPAbove1 2^Floor[exponent]];


(* ::Text:: *)
(*Extends the interval of its argument by a half ULP on both sides.  If Positive is True, the lower bound is not extended below 0.  This is pessimistic as we use the largest ULP over the interval.*)


Options[addHalfULPInterval]={Positive->False};
addHalfULPInterval[x_,OptionsPattern[]]:=
Block[{h=halfULP[x]},
(* The returned interval must explicitly contain x.  This ensures that, if x is unbound it can
later be removed to get the absolute error. *)
If[
OptionValue[Positive],
x+Interval[{-Min[h,IEEE754FloatingPointEvaluation`Private`intervalMin[x]],h}],
x+Interval[{-h,h}]]];


SetAttributes[IEEEEvaluateWithAbsoluteError,HoldAll];
Options[IEEEEvaluateWithAbsoluteError]={UseFMA->True};
IEEEEvaluateWithAbsoluteError[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evae,usefma=OptionValue[UseFMA]},
SetAttributes[evae,HoldAll];
evae[a_*b_+c_]:=If[
usefma,
addHalfULPInterval[evae[a]evae[b]+evae[c]],
addHalfULPInterval[addHalfULPInterval[evae[a]evae[b]]+evae[c]]];
evae[a_+b_]:=addHalfULPInterval[evae[a]+evae[b]];
evae[a_+b__]:=(Message[IEEEEvaluateWithAbsoluteError::badass]; $Failed);
evae[a_*b_]:=addHalfULPInterval[evae[a]evae[b]];
evae[a_*b__]:=(Message[IEEEEvaluateWithAbsoluteError::badass]; $Failed);
evae[a_/b_]:=addHalfULPInterval[evae[a]/evae[b]];
(* Negation is exact. *)
evae[-a_]:=-evae[a];
(* Squaring an interval is not the same as multiplying two identical intervals.
Also, if the lower bound of the square is 0, it is exact. *) 
evae[a_^2]:=addHalfULPInterval[evae[a]^2,Positive->True];
evae[a_^3]:=Block[{t},Expand[addHalfULPInterval[addHalfULPInterval[evae[t]^2,Positive->True]evae[t]]]/.t->a];
evae[a_^4]:=addHalfULPInterval[addHalfULPInterval[evae[a]^2,Positive->True]^2,Positive->True];
evae[a_?NumberQ]:=Block[{cra=CorrectlyRound[a]},Interval[{cra,cra}]];
evae[a_]:=ReleaseHold[a];
evae[x]];
IEEEEvaluateWithAbsoluteError[_, args__]:=
(Message[IEEEEvaluate::argnum, Length[{args}] + 1]; $Failed);


End[]


EndPackage[]
