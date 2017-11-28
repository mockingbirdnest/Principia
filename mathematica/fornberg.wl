(* ::Package:: *)

(* ::Title:: *)
(*Generation of Finite Difference Formulas on Arbitrary Spaced Grids*)


(* ::Subtitle:: *)
(*Bengt Fornberg*)


(* ::Subsubtitle:: *)
(*Mathematics of Computation, Volume 51, Number 184, doi:10.1090/S0025-5718-1988-0935077-0*)


(* ::Text:: *)
(*Implementation of the algorithm in section 2.*)


(* ::Input:: *)
(*GenerateFornberg[M_, N_, x0_, \[Alpha]_,OptionsPattern[]]:=*)
(*Module[*)
(*{\[Delta],c1=1},*)
(*\[Delta][0,0][0]=1;*)
(*Do[*)
(*Module[*)
(*{c2=1},*)
(*Do[*)
(*Module[*)
(*{c3=\[Alpha][n]-\[Alpha][\[Nu]]},*)
(*c2=c2 c3;*)
(*If[n<=M,\[Delta][n-1,\[Nu]][n]=0];*)
(*Do[*)
(*\[Delta][n,\[Nu]][m]=((\[Alpha][n]-x0)\[Delta][n-1,\[Nu]][m]-m \[Delta][n-1,\[Nu]][m-1])/c3,*)
(*{m,0,Min[n,M]}]*)
(*],*)
(*{\[Nu],0,n-1}*)
(*];*)
(*Do[*)
(*\[Delta][n,n][m]=(c1/c2)(m \[Delta][n-1,n-1][m-1]-(\[Alpha][n-1]-x0)\[Delta][n-1,n-1][m]),*)
(*{m,0,Min[n,M]}*)
(*];*)
(*c1=c2*)
(*],*)
(*{n,1,N}*)
(*];*)
(*Table[*)
(*Table[*)
(*If[*)
(*i<m-1,*)
(*0,*)
(*If[OddQ[m]&&OptionValue["Backwards"],-\[Delta][i,j][m],\[Delta][i,j][m]]*)
(*],*)
(*{i,0,N},*)
(*{j,0,i}*)
(*],*)
(*{m,0,M}*)
(*]*)
(*]*)


(* ::Input:: *)
(*Options[GenerateFornberg]={"Backwards"->False}*)


(* ::Text:: *)
(*A convenience function to reduce all the fractions to a common denominator and output the numerators and the denominator separately.*)


(* ::Input:: *)
(*GenerateSplitFornberg[M_, N_, x0_, \[Alpha]_,options___]:=*)
(*Module[*)
(*{*)
(*fornberg=GenerateFornberg[M,N,x0,\[Alpha],options],*)
(*commondenominators,*)
(*denominators*)
(*},*)
(*denominators=Map[Denominator[#]&,fornberg];*)
(*commondenominators=Map[Fold[LCM,#]&,denominators,{2}];*)
(*MapThread[{#1 #2,#2}&,{fornberg,commondenominators},2]*)
(*]*)


(* ::Input:: *)
(*GenerateFornberg[3,6,0,#&,"Backwards"->True]//TableForm*)


(* ::Input:: *)
(*GenerateSplitFornberg[1,6,0,#&,"Backwards"->True]//TableForm*)


(* ::Text:: *)
(*The data that we need for our integrators.*)


(* ::Input:: *)
(*GenerateSplitFornberg[1,14,0,#&,"Backwards"->True]*)
