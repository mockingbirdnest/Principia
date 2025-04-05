(* ::Package:: *)

(* ::Input:: *)
(*MullerWorstCaseRR[\[Beta]_,p_,eFirst_,eFinal_,C_,ndigits_]:=Module[{\[Epsilon]Min,powerOf\[Beta]OverC,e,a,pLast,r,qLast,Q,P,newQ,newP,\[Epsilon],numberMin,expMin,ell},*)
(*$MaxExtraPrecision=ndigits;*)
(*\[Epsilon]Min=12345;*)
(*powerOf\[Beta]OverC=\[Beta]^(eFirst-p)/C;*)
(*Do[*)
(*powerOf\[Beta]OverC=\[Beta] powerOf\[Beta]OverC;*)
(*a=Floor[powerOf\[Beta]OverC];*)
(*pLast=a;*)
(*r=1/(powerOf\[Beta]OverC-a);*)
(*a=Floor[r];*)
(*qLast=1;*)
(*Q=a;*)
(*P=pLast*a+1;*)
(*While[*)
(*Q<\[Beta]^p-1,*)
(*r=1/(r-a);*)
(*a=Floor[r];*)
(*newQ=Q a+qLast;*)
(*newP=P a+pLast;*)
(*qLast=Q;*)
(*pLast=P;*)
(*Q=newQ;*)
(*P=newP];*)
(*\[Epsilon]=N[C*Abs[pLast-qLast*powerOf\[Beta]OverC],ndigits];*)
(*If[\[Epsilon]<\[Epsilon]Min,\[Epsilon]Min=\[Epsilon];numberMin=qLast;expMin=e],*)
(*{e,eFirst-p+1,eFinal-p+1}];*)
(*Print[StringTemplate["The worst case occurs\n for x=``*``^``\n"][*)
(*numberMin,\[Beta],expMin]];*)
(*Print[StringTemplate["The corresponding reduced argument is:\n ``\n"][\[Epsilon]Min]];*)
(*ell=N[Log[\[Epsilon]Min]/Log[\[Beta]],10];*)
(*Print[StringTemplate["whose radix `` logarithm is ``"][\[Beta],ell]]];*)


(* ::Input:: *)
(*MullerWorstCaseRR[2,53,0,1023,\[Pi]/2,400]*)


(* ::Input:: *)
(*MullerWorstCaseRR[2,24,0,127,\[Pi]/2,400]*)


(* ::Input:: *)
(*MullerWorstCaseRR[2,27,15,19,\[Pi]/2,400]*)


(* ::Input:: *)
(*BaseForm[95255502,16]*)


(* ::Input:: *)
(*x=95255502*2^-8*)


(* ::Input:: *)
(*Log2[x]//N*)


(* ::Input:: *)
(*N[Sin[x],20]*)


(* ::Input:: *)
(*N[Cos[x],20]//Log2*)
