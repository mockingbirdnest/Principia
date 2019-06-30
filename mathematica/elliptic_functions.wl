(* ::Package:: *)

(* ::Section:: *)
(*Definitions*)


(* ::Input:: *)
(*decimalFloatLiteral[x_Real,exponentWidth_Integer]:=*)
(* With[*)
(*  {m=If[x==0,0,MantissaExponent[x][[1]]*10],*)
(*   e=If[x==0,0,MantissaExponent[x][[2]]-1]},*)
(*  StringJoin[*)
(*   StringRiffle[#,""]&/@*)
(*    {{#[[1]]},*)
(*     If[Length[#]>1,{"."},Nothing],*)
(*     If[Length[#]>1,StringPartition[#[[2]],UpTo[5]],Nothing],*)
(*     If[e!=0,"e"<>If[e>0,"+","-"]<>IntegerString[e,10,exponentWidth],Nothing]}&[*)
(*     StringSplit[ToString[m],"."]]]]*)


(* ::Input:: *)
(*decimalFloatLiteral[0,_]:="0"*)


(* ::Section:: *)
(*Random arguments with denser distribution near m=1*)


(* ::Input:: *)
(*SeedRandom[666]*)


(* ::Input:: *)
(*randomu=Sort[RandomReal[{-10,10},60,WorkingPrecision->20]];*)


(* ::Input:: *)
(*randomm=Sort[Select[1-RandomVariate[HalfNormalDistribution[3],60,WorkingPrecision->20],#>=0&]];*)


(* ::Input:: *)
(*randomargs=Flatten[Outer[List,randomu,randomm],1];*)


(* ::Input:: *)
(*randomvals=Map[*)
(*{*)
(*JacobiSN[#[[1]],#[[2]]],JacobiCN[#[[1]],#[[2]]],JacobiDN[#[[1]],#[[2]]],JacobiAmplitude[#[[1]],#[[2]]]}&,*)
(*randomargs];*)


(* ::Input:: *)
(*randomstrs=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[#[[1]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[2]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[3]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[4]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[5]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[6]],2]<>*)
(*"}"*)
(*&,*)
(*Join[randomargs,randomvals,2]];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Export[*)
(*"..\\numerics\\elliptic_functions.proto.txt",*)
(*StringRiffle[randomstrs,"\n"],*)
(*"text"]*)
