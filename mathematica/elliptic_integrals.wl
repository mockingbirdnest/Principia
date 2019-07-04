(* ::Package:: *)

(* ::Section:: *)
(*Definitions*)


(* ::Input:: *)
(*fukushimaB[\[CurlyPhi]_,m_]:=(EllipticE[\[CurlyPhi],m]-(1-m)EllipticF[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaD[\[CurlyPhi]_,m_]:=(EllipticF[\[CurlyPhi],m]-EllipticE[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaJ[\[CurlyPhi]_,n_,m_]:=(EllipticPi[n,\[CurlyPhi],m]-EllipticF[\[CurlyPhi],m])/n*)


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
(*Random arguments with denser distribution near \[CurlyPhi]=\[Pi]/2 and m=1*)


(* ::Subsection:: *)
(*Bivariate integrals*)


(* ::Input:: *)
(*SeedRandom[666]*)


(* ::Input:: *)
(*random\[CurlyPhi]=Sort[Select[\[Pi]/2-RandomVariate[HalfNormalDistribution[2],60,WorkingPrecision->20],#>=0&]];*)


(* ::Input:: *)
(*randomm=Sort[Select[1-RandomVariate[HalfNormalDistribution[3],60,WorkingPrecision->20],#>=0&]];*)


(* ::Input:: *)
(*randomargs2=Flatten[Outer[List,random\[CurlyPhi],randomm],1];*)


(* ::Input:: *)
(*randomvals2=Map[*)
(*{*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaB[#[[1]],m],m->0],*)
(*fukushimaB[#[[1]],#[[2]]]],*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaD[#[[1]],m],m->0],*)
(*fukushimaD[#[[1]],#[[2]]]],*)
(*EllipticE[#[[1]],#[[2]]],*)
(*EllipticF[#[[1]],#[[2]]]}&,*)
(*randomargs2];*)


(* ::Input:: *)
(*randomstrs2=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[#[[1]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[2]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[3]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[4]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[5]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[6]],2]<>*)
(*"}"*)
(*&,*)
(*Join[randomargs2,randomvals2,2]];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Export[*)
(*"..\\numerics\\bivariate_elliptic_integrals.proto.txt",*)
(*StringRiffle[randomstrs2,"\n"],*)
(*"text"]*)


(* ::Subsection:: *)
(*Trivariate integral*)


(* ::Input:: *)
(*SeedRandom[666]*)


(* ::Input:: *)
(*random\[CurlyPhi]=Sort[Select[\[Pi]/2-RandomVariate[HalfNormalDistribution[2],15,WorkingPrecision->20],#>=0&]];*)


(* ::Input:: *)
(*randomm=Sort[Select[1-RandomVariate[HalfNormalDistribution[3],15,WorkingPrecision->20],#>=0&]];*)


(* ::Input:: *)
(*randomn=Sort[RandomReal[{0,1},15,WorkingPrecision->20]];*)


(* ::Input:: *)
(*randomargs3=Flatten[Outer[List,random\[CurlyPhi],randomn,randomm],2];*)


(* ::Input:: *)
(*randomvals3=Map[*)
(*{*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaJ[#[[1]],n,#[[3]]],n->0],*)
(*fukushimaJ[#[[1]],#[[2]],#[[3]]]],*)
(*EllipticPi[#[[2]],#[[1]],#[[3]]]}&,*)
(*randomargs3];*)


(* ::Input:: *)
(*randomstrs3=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[#[[1]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[2]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[3]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[4]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[5]],2]<>*)
(*"}"*)
(*&,*)
(*Join[randomargs3,randomvals3,2]];*)


(* ::Subsubsection:: *)
(*A case that used to have a bug*)


(* ::Input:: *)
(*SeedRandom[666]*)


(* ::Input:: *)
(*random1\[CurlyPhi]=Sort[RandomReal[{0,\[Pi]/2},10,WorkingPrecision->20]];*)


(* ::Input:: *)
(*random1m=Sort[RandomReal[{999/1000,1},10,WorkingPrecision->20]];*)


(* ::Input:: *)
(*random1n=Sort[RandomReal[{0,1},10,WorkingPrecision->20]];*)


(* ::Input:: *)
(*random1args3=Flatten[Outer[List,random1\[CurlyPhi],random1n,random1m],2];*)


(* ::Input:: *)
(*random1vals3=Map[*)
(*{*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaJ[#[[1]],n,#[[3]]],n->0],*)
(*fukushimaJ[#[[1]],#[[2]],#[[3]]]],*)
(*EllipticPi[#[[2]],#[[1]],#[[3]]]}&,*)
(*random1args3];*)


(* ::Input:: *)
(*random1strs3=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[#[[1]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[2]],2]<>*)
(*" argument: "<>decimalFloatLiteral[#[[3]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[4]],2]<>*)
(*" value: "<>decimalFloatLiteral[#[[5]],2]<>*)
(*"}"*)
(*&,*)
(*Join[random1args3,random1vals3,2]];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Export[*)
(*"..\\numerics\\trivariate_elliptic_integrals.proto.txt",*)
(*StringRiffle[Flatten[{randomstrs3,random1strs3}],"\n"],*)
(*"text"]*)


(* ::Section:: *)
(*xeldbj near \[CurlyPhi]=\[Pi]/2*)


(* ::Input:: *)
(*xelbdj\[CurlyPhi]={14148475504056880/2^53}*)


(* ::Input:: *)
(*xelbdjn={9005001498122835/2^53,3/4,1/2,1/4,0}*)


(* ::Input:: *)
(*xelbdjm={1-8842084905851963/2^159,3/4,1/2,1/4,0}*)


(* ::Input:: *)
(*xelbdjargs=Flatten[Outer[List,xelbdjn,xelbdjm,xelbdj\[CurlyPhi]],2];*)


(* ::Input:: *)
(*xelbdjvals=Map[*)
(*N[{*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaB[#[[3]],m],m->0],*)
(*fukushimaB[#[[3]],#[[2]]]],*)
(*If[#[[2]]==0,*)
(*Limit[fukushimaD[#[[3]],m],m->0],*)
(*fukushimaD[#[[3]],#[[2]]]],*)
(*If[#[[1]]==0,*)
(*Limit[fukushimaJ[#[[3]],n,#[[2]]],n->0],*)
(*fukushimaJ[#[[3]],#[[1]],#[[2]]]]*)
(*},60]&,*)
(*xelbdjargs];*)


(* ::Input:: *)
(*xelbdjstrs=Map[*)
(*"entry { argument: "<>decimalFloatLiteral[N[#[[1]],5],2]<>*)
(*" argument: "<>decimalFloatLiteral[N[#[[2]],5],2]<>*)
(*" argument: "<>decimalFloatLiteral[N[#[[3]]/\[Pi],5],2]<>*)
(*" value: "<>decimalFloatLiteral[N[#[[4]],21],2]<>*)
(*" value: "<>decimalFloatLiteral[N[#[[5]],21],2]<>*)
(*" value: "<>decimalFloatLiteral[N[#[[6]],21],2] <>*)
(*"}"*)
(*&,*)
(*Join[xelbdjargs,xelbdjvals,2]];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Input:: *)
(*Export[*)
(*"..\\temp\\xelbdj_\[Pi]_over_2.proto.txt",*)
(*StringRiffle[xelbdjstrs,"\n"],*)
(*"text"]*)
