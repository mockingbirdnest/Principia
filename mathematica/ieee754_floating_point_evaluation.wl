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
"evaluates \!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\) using the rules of " <>
"IEEE arithmetic, and propagates the absolute error bounds at each stage of " <>
"the computation.  Returns a list of two intervals: the first one is the " <>
"interval of the result evaluated with proper IEEE rounding, the second one is " <>
"the range of the absolute error on the result.  The argument is an expression " <>
"which can include numbers (assumed to be exact after correct rounding), " <>
"intervals (assumed to not carry any error), or a list of two intervals " <> 
"for a value and its absolute error."
IEEEEvaluateWithAbsoluteError::argnum =
"IEEEEvaluateWithAbsoluteError called with `1` arguments; 1 argument is expected.";
IEEEEvaluateWithAbsoluteError::badass =
"IEEEEvaluateWithAbsoluteError does not support associativity, expressions must be " <>
"parenthesized."


IEEEEvaluateWithRelativeError;
IEEEEvaluateWithRelativeError::usage =
"IEEEEvaluateWithRelativeeError[\!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\)] " <>
"evaluates \!\(\*StyleBox[\"x\",\nFontSlant->\"Italic\"]\) using the rules of " <>
"IEEE arithmetic, and propagates the relative error bounds at each stage of " <>
"the computation.  Returns a list of two intervals: the first one is the " <>
"interval of the result evaluated with proper IEEE rounding, the second one is " <>
"the range of the relative error on the result.  The argument is an expression " <>
"which can include numbers (assumed to be exact after correct rounding), " <>
"intervals (assumed to not carry any error), or a list of two intervals " <>
"for a value and its relative error.";
IEEEEvaluateWithRelativeError::argnum =
"IEEEEvaluateWithRelativeError called with `1` arguments; 1 argument is expected.";
IEEEEvaluateWithRelativeError::badass =
"IEEEEvaluateWithRelativeError does not support associativity, expressions must be " <>
"parenthesized."


UseFMA;
UseFMA::usage =
"UseFMA is an option for IEEEEvaluate, IEEEEvaluateWithAbsoluteError, and " <>
"IEEEEvaluateWithRelativeError that specifies whether to use FMA for " <>
"expressions of the form \!\(\*
StyleBox[\"a\",\nFontSlant->\"Italic\"]\) * \!\(\*
StyleBox[\"b\",\nFontSlant->\"Italic\"]\) + \!\(\*
StyleBox[\"c\",\nFontSlant->\"Italic\"]\).";


Exact;
Exact::usage =
"Exact is a wrapper for expressions passed to IEEEEvaluateWithAbsoluteError, and " <>
"IEEEEvaluateWithRelativeError that specifies that the top-level operation " <>
"is exact (errors may still be propagated from the subexpressions."


Begin["`Private`"]


SetAttributes[Exact,HoldAll];


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


(* ::Text:: *)
(*The error bound on an IEEE computation in the binade [1/2, 1[.*)


errorBelow1:=If[
	RoundingMode[]==NearestTiesToEven,
	1-FromRepresentation[Representation[1]-1/2],
	1-FromRepresentation[Representation[1]-1]];


(* ::Text:: *)
(*Returns the error bound for the largest element (in absolute value) of its argument.  The returned value is positive, regardless of the sign of the argument.  The argument is an interval.  Note that if the largest element is a power of two, the error bound is the one below that power of two.*)


(* Would want to write Max[Log2[...]] below but that doesn't work somehow. *)
absoluteErrorBound[x_]:=Block[{exponent=Log2[Max[Abs[x]]]},errorBelow1 2^Ceiling[exponent]];


Options[applyOpWithAbsoluteError]={Exact->False};
(* Special case because for an interval x*x^2 is not x^3. *)
applyOpWithAbsoluteError[cube,{va_,\[Delta]a_},OptionsPattern[]]:=Block[
	{va2,\[Delta]a2,h,r,\[Delta]r},
	{va2,\[Delta]a2}=applyOpWithAbsoluteError[#^2&,{va,\[Delta]a}];
	r=va^3;(* Wrong! *)
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,absoluteErrorBound[r]];
	\[Delta]r=ReplaceAll[
		Expand[(v2+\[Delta]2)(v+\[Delta])-(v2 v)],
		{v->va,\[Delta]->\[Delta]a,v2->va2,\[Delta]2->\[Delta]a2}]+Interval[{-h,h}];
	{r,\[Delta]r}];
applyOpWithAbsoluteError[op_,{va_,\[Delta]a_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,absoluteErrorBound[r]];
	\[Delta]r=ReplaceAll[Expand[op[v+\[Delta]]-op[v]],{v->va,\[Delta]->\[Delta]a}]+Interval[{-h,h}];
	{r,\[Delta]r}];
applyOpWithAbsoluteError[op_,{va_,\[Delta]a_},{vb_,\[Delta]b_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va,vb];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,absoluteErrorBound[r]];
	\[Delta]r=ReplaceAll[
		Expand[op[v1+\[Delta]1,v2+\[Delta]2]-op[v1,v2]],
		{v1->va,\[Delta]1->\[Delta]a,v2->vb,\[Delta]2->\[Delta]b}]+Interval[{-h,h}];
	{r,\[Delta]r}];
applyOpWithAbsoluteError[op_,{va_,\[Delta]a_},{vb_,\[Delta]b_},{vc_,\[Delta]c_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va,vb,vc];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,absoluteErrorBound[r]];
	\[Delta]r=ReplaceAll[
		Expand[op[v1+\[Delta]1,v2+\[Delta]2,v3+\[Delta]3]-op[v1,v2,v3]],
		{v1->va,\[Delta]1->\[Delta]a,v2->vb,\[Delta]2->\[Delta]b,v3->vc,\[Delta]3->\[Delta]c}]+Interval[{-h,h}];
	{r,\[Delta]r}];


SetAttributes[IEEEEvaluateWithAbsoluteError,HoldAll];
Options[IEEEEvaluateWithAbsoluteError]={UseFMA->True};
IEEEEvaluateWithAbsoluteError[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evae,usefma=OptionValue[UseFMA]},
SetAttributes[evae,HoldAll];
Options[evae]={Exact->False};
evae[a_*b_+c_,opts:OptionsPattern[]]:=If[
usefma,
applyOpWithAbsoluteError[#1 #2+#3&,evae[a],evae[b],evae[c],opts],
applyOpWithAbsoluteError[
	Plus,applyOpWithAbsoluteError[Times,evae[a],evae[b],opts],evae[c]],opts];
evae[a_+b_,opts:OptionsPattern[]]:=applyOpWithAbsoluteError[Plus,evae[a],evae[b],opts];
evae[a_+b__]:=(Message[IEEEEvaluateWithAbsoluteError::badass]; $Failed);
evae[a_*b_,opts:OptionsPattern[]]:=applyOpWithAbsoluteError[Times,evae[a],evae[b],opts];
evae[a_*b__]:=(Message[IEEEEvaluateWithAbsoluteError::badass]; $Failed);
evae[a_/b_,opts:OptionsPattern[]]:=applyOpWithAbsoluteError[Divide,evae[a],evae[b],opts];
(* Negation is exact. *)
evae[-a_]:=-evae[a];
(* Squaring an interval is not the same as multiplying two identical intervals. *) 
evae[a_^2,opts:OptionsPattern[]]:=applyOpWithAbsoluteError[#^2&,evae[a],opts];
evae[a_^3,opts:OptionsPattern[]]:=applyOpWithAbsoluteError[cube,evae[a],opts];
evae[a_^4,opts:OptionsPattern[]]:=
	applyOpWithAbsoluteError[#^2&,applyOpWithAbsoluteError[#^2&,evae[a],opts],opts];
evae[a_?NumberQ]:=Block[{cra=CorrectlyRound[a]},evae[Interval[{cra,cra}]]];
evae[{v_Interval,\[Delta]_Interval}]:={v,\[Delta]};
evae[a_Interval]:={a,Interval[{0,0}]};
evae[Exact[a_]]:=evae[a,Exact->True];
evae[a_?ValueQ]:=evae[Evaluate[a]];
evae[a_]:=a;
evae[x]];
IEEEEvaluateWithAbsoluteError[_, args__]:=
(Message[IEEEEvaluateWithAbsoluteError::argnum, Length[{args}] + 1]; $Failed);


(* ::Text:: *)
(*The relative error bound on an IEEE computation.*)


relativeErrorBound := If[
	RoundingMode[]==NearestTiesToEven,
	FromRepresentation[Representation[1]+1/2]-1,
	FromRepresentation[Representation[1]+1]-1];


(* ::Text:: *)
(*Returns the relative error bound for a computation.  The relative error is "small", i.e., close to 2^-53, not to (1+2^-53), and is an interval.  The arguments are pairs of {value interval, relative error interval}.*)


Options[applyOpWithRelativeError]={Exact->False};
(* Special case because for an interval x*x^2 is not x^3. *)
applyOpWithRelativeError[cube,{va_,\[Delta]a_},OptionsPattern[]]:=Block[
	{va2,\[Delta]a2,h,r,\[Delta]r},
	{va2,\[Delta]a2}=applyOpWithRelativeError[#^2&,{va,\[Delta]a}];
	r=va^3;
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,relativeErrorBound];
	\[Delta]r=Interval[{1-h,1+h}]
		(ReplaceAll[Expand[Apart[(v2(1+\[Delta]2)v(1+\[Delta]))/(v2 v)-1]],
			{v->va,\[Delta]->\[Delta]a,v2->va2,\[Delta]2->\[Delta]a2}]+1)-1;
	{r,\[Delta]r}];
applyOpWithRelativeError[op_,{va_,\[Delta]a_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,relativeErrorBound];
	\[Delta]r=Interval[{1-h,1+h}]
		(ReplaceAll[Expand[Apart[op[v(1+\[Delta])]/op[v]-1]],{v->va,\[Delta]->\[Delta]a}]+1)-1;
	{r,\[Delta]r}];
applyOpWithRelativeError[op_,{va_,\[Delta]a_},{vb_,\[Delta]b_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va,vb];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,relativeErrorBound];
	\[Delta]r=Interval[{1-h,1+h}]
		(ReplaceAll[Expand[Apart[op[v1(1+\[Delta]1),v2(1+\[Delta]2)]/op[v1,v2]-1]],
			{v1->va,\[Delta]1->\[Delta]a,v2->vb,\[Delta]2->\[Delta]b}]+1)-1;
	{r,\[Delta]r}];
applyOpWithRelativeError[op_,{va_,\[Delta]a_},{vb_,\[Delta]b_},{vc_,\[Delta]c_},OptionsPattern[]]:=Block[
	{h,r,\[Delta]r},
	r=op[va,vb,vc];
	r=Interval[{CorrectlyRound[Min[r]],CorrectlyRound[Max[r]]}];
	h=If[OptionValue[Exact],0,relativeErrorBound];
	\[Delta]r=Interval[{1-h,1+h}]
		(ReplaceAll[Expand[Apart[op[v1(1+\[Delta]1),v2(1+\[Delta]2),v3(1+\[Delta]3)]/op[v1,v2,v3]-1]],
			{v1->va,\[Delta]1->\[Delta]a,v2->vb,\[Delta]2->\[Delta]b,v3->vc,\[Delta]3->\[Delta]c}]+1)-1;
	{r,\[Delta]r}];


SetAttributes[IEEEEvaluateWithRelativeError,HoldAll];
Options[IEEEEvaluateWithRelativeError]={UseFMA->True};
IEEEEvaluateWithRelativeError[x_,OptionsPattern[]]:=
Block[
{Plus,Times,evre,usefma=OptionValue[UseFMA]},
SetAttributes[evre,HoldAll];
Options[evre]={Exact->False};
evre[a_*b_+c_,opts:OptionsPattern[]]:=If[
usefma,
applyOpWithRelativeError[
	#1 #2+#3&,evre[a],evre[b],evre[c],opts],
applyOpWithRelativeError[
	Plus,applyOpWithRelativeError[Times,evre[a],evre[b],opts],evre[c]],opts];
evre[a_+b_,opts:OptionsPattern[]]:=applyOpWithRelativeError[Plus,evre[a],evre[b],opts];
evre[a_+b__]:=(Message[IEEEEvaluateWithRelativeError::badass]; $Failed);
evre[a_*b_,opts:OptionsPattern[]]:=applyOpWithRelativeError[Times,evre[a],evre[b],opts];
evre[a_*b__]:=(Message[IEEEEvaluateWithRelativeError::badass]; $Failed);
evre[a_/b_,opts:OptionsPattern[]]:=applyOpWithRelativeError[Divide,evre[a],evre[b],opts];
(* Negation is exact. *)
evre[-a_]:=-evre[a];
(* Squaring an interval is not the same as multiplying two identical intervals. *) 
evre[a_^2,opts:OptionsPattern[]]:=applyOpWithRelativeError[#^2&,evre[a],opts];
evre[a_^3,opts:OptionsPattern[]]:=applyOpWithRelativeError[cube,evre[a],opts];
evre[a_^4,opts:OptionsPattern[]]:=
	applyOpWithRelativeError[#^2&,applyOpWithRelativeError[#^2&,evre[a],opts],opts];
evre[a_?NumberQ]:=Block[{cra=CorrectlyRound[a]},evre[Interval[{cra,cra}]]];
evre[{v_Interval,\[Delta]_Interval}]:={v,\[Delta]};
evre[a_Interval]:={a,Interval[{0,0}]};
evre[Exact[a_]]:=evre[a,Exact->True];
evre[a_?ValueQ]:=evre[Evaluate[a]];
evre[a_]:=a;
evre[x]];
IEEEEvaluateWithRelativeError[_, args__]:=
(Message[IEEEEvaluateWithRelativeError::argnum, Length[{args}] + 1]; $Failed);


End[]


EndPackage[]
