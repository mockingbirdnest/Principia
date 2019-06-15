(* ::Package:: *)

(* ::Input:: *)
(*fukushimaB[\[CurlyPhi]_,m_]:=(EllipticE[\[CurlyPhi],m]-(1-m)EllipticF[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaD[\[CurlyPhi]_,m_]:=(EllipticF[\[CurlyPhi],m]-EllipticE[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaJ[\[CurlyPhi]_,n_,m_]:=(EllipticPi[n,\[CurlyPhi],m]-EllipticF[\[CurlyPhi],m])/n*)


(* ::Input:: *)
(*decimalFloatLiteral[x_Real,exponentWidth_Integer]:=*)
(* With[*)
(*  {m=MantissaExponent[x][[1]]*10,*)
(*   e=MantissaExponent[x][[2]]-1},*)
(*  StringJoin[*)
(*   StringRiffle[#,""]&/@*)
(*    {{#[[1]]},*)
(*     If[Length[#]>1,{"."},Nothing],*)
(*     If[Length[#]>1,StringPartition[#[[2]],UpTo[5]],Nothing],*)
(*     If[e!=0,"e"<>If[e>0,"+","-"]<>IntegerString[e,10,exponentWidth],Nothing]}&[*)
(*     StringSplit[ToString[m],"."]]]]*)


(* ::Input:: *)
(*SeedRandom[666]*)


(* ::Input:: *)
(*random\[CurlyPhi]=Sort[RandomReal[{0,\[Pi]/2},10,WorkingPrecision->20]]*)


(* ::Input:: *)
(*randomn=Sort[RandomReal[{0,1},10,WorkingPrecision->20]]*)


(* ::Input:: *)
(*randomm=Sort[RandomReal[{999/1000,1},10,WorkingPrecision->20]]*)


(* ::Input:: *)
(*args=Flatten[Outer[List,randomn,randomm,random\[CurlyPhi]],2];*)


(* ::Input:: *)
(*vals=Map[*)
(*{*)
(*#[[1]],#[[2]],#[[3]],*)
(*fukushimaB[#[[3]],#[[2]]],*)
(*fukushimaD[#[[3]],#[[2]]],*)
(*fukushimaJ[#[[3]],#[[1]],#[[2]]]*)
(*}&,*)
(*args,*)
(*{1}];*)


(* ::Input:: *)
(*strs=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[#[[1]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[2]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[3]],2]<>*)
(*" value: "<>decimalFloatLiteral[fukushimaB[#[[3]],#[[2]]],2]<>*)
(*" value: "<>decimalFloatLiteral[fukushimaD[#[[3]],#[[2]]],2]<>*)
(*" value: "<>decimalFloatLiteral[fukushimaJ[#[[3]],#[[1]],#[[2]]],2] <>*)
(*"}"*)
(*&,*)
(*args,*)
(*{1}];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Export[*)
(*"..\\numerics\\elliptic_integrals.proto.txt",*)
(*StringRiffle[strs,"\n"],*)
(*"text"]*)
