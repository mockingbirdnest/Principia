(* ::Package:: *)

(* ::Title:: *)
(*Elements of the Outer Planets for One Million Years*)


(* ::Subtitle:: *)
(*C. J. Cohen, E. C. Hubbard, Claus Oesterwinter*)


(* ::Subsubtitle:: *)
(*Astronomical papers prepared for the use of the American ephemeris and nautical almanac, Volume XXII, Part I*)


SetDirectory[NotebookDirectory[]];
<<"fornberg.wl"


fs[x_,n_]:=Normal[Series[f[\[Xi]],{\[Xi],0,n}]]/.\[Xi]->x;
f2s[x_,n_]:=Normal[Series[f''[\[Xi]],{\[Xi],0,n}]]/.\[Xi]->x;


(* ::Text:: *)
(*The \[Beta] coefficients from page 21.*)


\[Beta][\[ScriptCapitalN]_]:=Module[
{\[CapitalDelta] = GenerateFornberg[\[ScriptCapitalN], \[ScriptCapitalN], 1, -# &],\[Delta]},
\[Delta][m_, n_, j_] := \[CapitalDelta][[m + 1, n + 1, j + 1]];
With[
 {n = \[ScriptCapitalN]},
 Table[
  Sum[
   ((-1)^(k - 1) (1 - 2^k) \[Delta][k - 2, n - 2, j])/k!,
   {k, 2, n}],
  {j, 0, n - 2}]]];


{{LCM@@(Denominator/@#),#*LCM@@(Denominator/@#)}}&[\[Beta][14]]//TableForm


With[{\[ScriptCapitalN]=14},(fs[-h,\[ScriptCapitalN]+3]-fs[-2h,\[ScriptCapitalN]+3])/h+h Sum[f2s[-i h,\[ScriptCapitalN]+3]\[Beta][\[ScriptCapitalN]][[i]],{i,1,Length[\[Beta][\[ScriptCapitalN]]]}]]//Expand


(* ::Text:: *)
(*A similar formula, using f(0) and f(-h) instead of f(-h) and f(-2h).*)


\[Gamma][\[ScriptCapitalN]_]:=Module[
{\[CapitalDelta] = GenerateFornberg[\[ScriptCapitalN], \[ScriptCapitalN], 1, -# &],\[Delta]},
\[Delta][m_, n_, j_] := \[CapitalDelta][[m + 1, n + 1, j + 1]];
With[
 {n = \[ScriptCapitalN]},
 Table[
  Sum[
   (-1)^k(\[Delta][k - 2, n - 2, j])/k!,
   {k, 2, n}],
  {j, 0, n - 2}]]];


With[{\[ScriptCapitalN]=14},(fs[0,\[ScriptCapitalN]+3]-fs[-h,\[ScriptCapitalN]+3])/h+h Sum[f2s[-i h,\[ScriptCapitalN]+3]\[Gamma][\[ScriptCapitalN]][[i]],{i,1,Length[\[Gamma][\[ScriptCapitalN]]]}]]//Expand


(* ::Text:: *)
(*A similar formula, using f(0) and f(-h) instead of f(-h) and f(-2h), and f''(0) through f''(-12h) instead of f''(-h) through f''(-13h).*)


\[Eta][\[ScriptCapitalN]_]:=Module[
{\[CapitalDelta] = GenerateFornberg[\[ScriptCapitalN], \[ScriptCapitalN], 0, -# &],\[Delta]},
\[Delta][m_, n_, j_] := \[CapitalDelta][[m + 1, n + 1, j + 1]];
With[
 {n = \[ScriptCapitalN]},
 Table[
  Sum[
   (-1)^k(\[Delta][k - 2, n - 2, j])/k!,
   {k, 2, n}],
  {j, 0, n - 2}]]];


With[{\[ScriptCapitalN]=14},(fs[0,\[ScriptCapitalN]+3]-fs[-h,\[ScriptCapitalN]+3])/h+h Sum[f2s[-(i-1) h,\[ScriptCapitalN]+3]\[Eta][\[ScriptCapitalN]][[i]],{i,1,Length[\[Eta][\[ScriptCapitalN]]]}]]//Expand


(* ::Text:: *)
(*The backward difference formula.*)


\[Delta][\[ScriptCapitalN]_]:=GenerateFornberg[\[ScriptCapitalN], \[ScriptCapitalN], 0, -# &][[1+1,\[ScriptCapitalN]+1]];


With[{\[ScriptCapitalN]=14},1/h Sum[fs[-(i-1) h,\[ScriptCapitalN]+3]\[Delta][\[ScriptCapitalN]][[i]],{i,1,Length[\[Delta][\[ScriptCapitalN]]]}]]//Expand


(* ::Text:: *)
(*The central difference formula; note that this does not actually use f(0).*)


\[Kappa][\[ScriptCapitalN]_]:=GenerateFornberg[\[ScriptCapitalN], \[ScriptCapitalN], 0, (-1)^# 2Ceiling[#/2] &][[1+1,\[ScriptCapitalN]+1]];


With[{\[ScriptCapitalN]=14},1/h Sum[fs[(-1)^(i-1)2Ceiling[(i-1)/2] h,\[ScriptCapitalN]+3]\[Kappa][\[ScriptCapitalN]][[i]],{i,1,Length[\[Kappa][\[ScriptCapitalN]]]}]]//Expand


(* ::Text:: *)
(*Comparing the coefficient of the lowest-order error term.*)


293074482523 /1120863744000//N


1611122953987 /7846046208000//N


45183033541/15692092416000//N


1/15//N


2048 /6435//N


GenerateCohenHubbardOesterwinter[\[Alpha]_, \[ScriptCapitalN]_]:=
Module[
{cho=Table[\[Alpha][n],{n,2,\[ScriptCapitalN]}],
commondenominators,
denominators},
denominators=Map[Denominator[#]&,cho];
commondenominators=Map[Fold[LCM,#]&,denominators,{1}];
MapThread[{#1 #2,#2}&,{cho,commondenominators},1]
]


g=GenerateCohenHubbardOesterwinter[\[Eta],14]
