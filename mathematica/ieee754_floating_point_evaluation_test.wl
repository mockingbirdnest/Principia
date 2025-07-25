(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


(* ::Input::Initialization:: *)
Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]];Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point_evaluation.wl"}]];


(* ::Input::Initialization:: *)
SetFloatingPointFormat[binary64]
SetRoundingMode[NearestTiesToEven]


(* ::Input::Initialization:: *)
On[Assert]


(* ::Input::Initialization:: *)
?"IEEE754FloatingPointEvaluation`*"


(* ::Section:: *)
(*IEEEEvaluate*)


(* ::Input::Initialization:: *)
?UseFMA


(* ::Input::Initialization:: *)
?IEEEEvaluate


(* ::Input::Initialization:: *)
Assert[RoundingMode[]==NearestTiesToEven];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1]==1];
Assert[IEEEEvaluate[1+2]==3];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1.5]==3/2];
Assert[IEEEEvaluate[0.1`100]==CorrectlyRound[1/10]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3]==CorrectlyRound[1/3]];
Assert[IEEEEvaluate[1/3+1/3]==2CorrectlyRound[1/3]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3+1/5]==CorrectlyRound[CorrectlyRound[1/3]+CorrectlyRound[1/5]]];
Assert[IEEEEvaluate[1/3-1/5]==CorrectlyRound[CorrectlyRound[1/3]-CorrectlyRound[1/5]]];
Assert[IEEEEvaluate[(1/3)*(1/5)]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[(1/3)^2]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]];
Assert[IEEEEvaluate[(1/3)^3]==
CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[1/3]]];
Assert[IEEEEvaluate[(1/3)^4]==
CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]*CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/3]]]];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[N[1/3,100]]==CorrectlyRound[1/3]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[\[Pi]]==$Failed]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]==CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]+CorrectlyRound[1/7]]];
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]==CorrectlyRound[CorrectlyRound[CorrectlyRound[1/3]*CorrectlyRound[1/5]]+CorrectlyRound[1/7]]];
Assert[IEEEEvaluate[1/3*1/5+1/7,UseFMA->False]!=IEEEEvaluate[1/3*1/5+1/7,UseFMA->True]];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((1/10+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10]==
Nest[CorrectlyRound[#+CorrectlyRound[1/10]]&,0,10]];
Assert[
IEEEEvaluate[((((((((1/10+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10)+1/10]!=1];
Assert[
IEEEEvaluate[1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10+1/10]==$Failed];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((1/10*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10]==
Nest[CorrectlyRound[#*CorrectlyRound[1/10]]&,1,10]];
Assert[
IEEEEvaluate[((((((((1/10*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10)*1/10]!=10^-10];
Assert[
IEEEEvaluate[1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10*1/10]==$Failed];


(* ::Input::Initialization:: *)
Assert[
IEEEEvaluate[((((((((-1/10-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10]==
Nest[CorrectlyRound[#-CorrectlyRound[1/10]]&,0,10]];
Assert[
IEEEEvaluate[((((((((-1/10-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10)-1/10]!=-1];
Assert[
IEEEEvaluate[-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10-1/10]==$Failed];


(* ::Input::Initialization:: *)
Module[
{x=3},
Assert[IEEEEvaluate[x]==3]];


(* ::Input::Initialization:: *)
Module[
{a,b,poly,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly[u_]:=u^3(a+b u^2);
Assert[IEEEEvaluate[poly[v]]!=CorrectlyRound[poly[v]]];
Assert[IEEEEvaluate[poly[v]]==
CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
CorrectlyRound[
CorrectlyRound[a]+
CorrectlyRound[b]CorrectlyRound[w^2]]]]]


(* ::Input::Initialization:: *)
Module[
{a,b,poly,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly=Function[u,u^3(a+b u^2)];
Assert[IEEEEvaluate[poly[v]]!=CorrectlyRound[poly[v]]];
Assert[IEEEEvaluate[poly[v]]==
CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
CorrectlyRound[
CorrectlyRound[a]+
CorrectlyRound[b]CorrectlyRound[w^2]]]]]


(* ::Section:: *)
(*IEEEEvaluateWithAbsoluteError*)


(* ::Input::Initialization:: *)
?IEEEEvaluateWithAbsoluteError


(* ::Input::Initialization:: *)
$ContextPath=Prepend[$ContextPath,"IEEE754FloatingPointEvaluation`Private`"]


(* ::Input::Initialization:: *)
Assert[errorBelow1==2^-54];


(* ::Input::Initialization:: *)
Assert[absoluteErrorBound[Interval[{1.2,1.3}]]==2^-53];
Assert[absoluteErrorBound[Interval[{-2.4,1.3}]]==2^-52];
Assert[absoluteErrorBound[Interval[{1,1}]]==2^-54];
Assert[absoluteErrorBound[Interval[{-1/2,1/2}]]==2^-55];


(* ::Input::Initialization:: *)
$ContextPath=Drop[$ContextPath,1]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[2*3+4,UseFMA->True]=={Interval[{10,10}],Interval[{-2^-50,2^-50}]}];
Assert[IEEEEvaluateWithAbsoluteError[2*3+4,UseFMA->False]=={Interval[{10,10}],Interval[{-2^-51-2^-50,2^-51+2^-50}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[1+2]=={Interval[{3,3}],Interval[{-2^-52,2^-52}]}];
Assert[IEEEEvaluateWithAbsoluteError[2*3]=={Interval[{6,6}],Interval[{-2^-51,2^-51}]}];
Assert[IEEEEvaluateWithAbsoluteError[2/3]=={Interval[{CorrectlyRound[2/3],CorrectlyRound[2/3]}],Interval[{-2^-54,2^-54}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[-Interval[{-2,1}]]=={Interval[{-1,2}],Interval[{0,0}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,2}]^2]=={Interval[{0,4}],Interval[{-2^-52,2^-52}]}];
(* The error after the squaring is [0, 2^-50]. *)
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^3]=={Interval[{-1,27}],Interval[{-3 2^-50-2^-49,3 2^-50+2^-49}]}];
Assert[IEEEEvaluateWithAbsoluteError[Interval[{-1,3}]^4]=={Interval[{0,81}],Interval[{-18 2^-50-2^-47,2^-100+18 2^-50+2^-47}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[CorrectlyRound[\[Pi]]]=={Interval[{884279719003555/281474976710656,884279719003555/281474976710656}],Interval[{0,0}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[
	{Interval[{2,2}],Interval[{-2^-50,2^-50}]}+
	{Interval[{3,3}],Interval[{-2^-51,2^-51}]}]==
	{Interval[{5,5}],Interval[{-2^-50-2 2^-51,2^-50+2 2^-51}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[
Exact[
		{Interval[{2,2}],Interval[{-2^-50,2^-50}]}+
		{Interval[{3,3}],Interval[{-2^-51,2^-51}]}]]==
	{Interval[{5,5}],Interval[{-2^-50- 2^-51,2^-50+ 2^-51}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[\[Pi]]==$Failed];
Assert[IEEEEvaluateWithAbsoluteError[33333333333333333]==$Failed];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[0.1`20]=={Interval[{3602879701896397/36028797018963968,3602879701896397/36028797018963968}],Interval[{3602879701896397/36028797018963968-0.1`20,3602879701896397/36028797018963968-0.1`20}]}]


(* ::Input::Initialization:: *)
Module[
{x=3},
Assert[IEEEEvaluateWithAbsoluteError[x]=={Interval[{3,3}],Interval[{0,0}]}]];


(* ::Input::Initialization:: *)
Module[
{a,b,poly,p,q,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly[u_]:=u^3(a+b u^2);
p=CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
	CorrectlyRound[
	CorrectlyRound[a]+
	CorrectlyRound[b]CorrectlyRound[w^2]]];
q=CorrectlyRound[poly[v]];
Assert[IEEEEvaluateWithAbsoluteError[poly[v]][[1]]!=Interval[{q,q}]];
Assert[IEEEEvaluateWithAbsoluteError[poly[v]][[1]]==Interval[{p,p}]]]


(* ::Input::Initialization:: *)
Module[
{a,b,poly,p,q,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly=Function[u,u^3(a+b u^2)];
p=CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
	CorrectlyRound[
	CorrectlyRound[a]+
	CorrectlyRound[b]CorrectlyRound[w^2]]];
q=CorrectlyRound[poly[v]];
Assert[IEEEEvaluateWithAbsoluteError[poly[v]][[1]]!=Interval[{q,q}]];
Assert[IEEEEvaluateWithAbsoluteError[poly[v]][[1]]==Interval[{p,p}]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithAbsoluteError[2,3]==$Failed];


(* ::Subsection:: *)
(*[SZ05] Analysis*)


(* ::Text:: *)
(*The error analysis from [SZ05], section 2.1 (section 3.1 in the preprint, which is more explicit when it comes to the analysis.)*)


(* ::Input::Initialization:: *)
Begin["`SZ05`"]


(* ::Input::Initialization:: *)
bits[{v_Interval,\[Delta]v_Interval}]:=N[Log2[Map[Max,Abs[{v,\[Delta]v}]]]];


(* ::Input::Initialization:: *)
a2=1/2


(* ::Input::Initialization:: *)
a3=CorrectlyRound[1/3!]


(* ::Input::Initialization:: *)
a4=CorrectlyRound[1/4!]


(* ::Input::Initialization:: *)
a5=CorrectlyRound[1/5!]


(* ::Input::Initialization:: *)
h=Interval[{CorrectlyRound[-(2^-10.977``100)],CorrectlyRound[2^-10.977``100]}]


(* ::Input::Initialization:: *)
k=IEEEEvaluateWithAbsoluteError[h^2]


(* ::Input::Initialization:: *)
bits[k]


(* ::Input::Initialization:: *)
k\[Prime]=IEEEEvaluateWithAbsoluteError[h^3]


(* ::Input::Initialization:: *)
bits[k\[Prime]]


(* ::Input::Initialization:: *)
S\[FivePointedStar]1=IEEEEvaluateWithAbsoluteError[a5 k]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]1]


(* ::Input::Initialization:: *)
S\[FivePointedStar]2=IEEEEvaluateWithAbsoluteError[S\[FivePointedStar]1-a3]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]2]


(* ::Input::Initialization:: *)
S\[FivePointedStar]3=IEEEEvaluateWithAbsoluteError[k\[Prime] S\[FivePointedStar]2]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]3]


(* ::Input::Initialization:: *)
C\[FivePointedStar]1=IEEEEvaluateWithAbsoluteError[a4 k]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]1]


(* ::Input::Initialization:: *)
C\[FivePointedStar]2=IEEEEvaluateWithAbsoluteError[C\[FivePointedStar]1-a2]


(* ::Text:: *)
(*NOTE: [SZ05] has an error on the 56th bit below, but their calculation looks wrong because the interval doesn't have bounds below 1/4.*)


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]2]


(* ::Input::Initialization:: *)
C\[FivePointedStar]3=IEEEEvaluateWithAbsoluteError[k C\[FivePointedStar]2]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]3]


(* ::Text:: *)
(*In a single operation :*)


(* ::Input::Initialization:: *)
S\[FivePointedStar]=IEEEEvaluateWithAbsoluteError[(a5 h^2-a3)h^3]


(* ::Input::Initialization:: *)
bits[S\[FivePointedStar]]


(* ::Input::Initialization:: *)
C\[FivePointedStar]=IEEEEvaluateWithAbsoluteError[h^2(a4 h^2-a2)]


(* ::Input::Initialization:: *)
bits[C\[FivePointedStar]]


(* ::Input::Initialization:: *)
End[]


(* ::Section:: *)
(*IEEEEvaluateWithRelativeError*)


?IEEEEvaluateWithRelativeError


(* ::Input::Initialization:: *)
$ContextPath=Prepend[$ContextPath,"IEEE754FloatingPointEvaluation`Private`"]


(* ::Input::Initialization:: *)
Assert[relativeErrorBound==2^-53];


(* ::Input::Initialization:: *)
$ContextPath=Drop[$ContextPath,1]


(* ::Subsection:: *)
(*Explicit Bounds*)


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[2*3+4,UseFMA->True]=={Interval[{10,10}],Interval[{-2^-53,2^-53}]}];
Assert[IEEEEvaluateWithRelativeError[2*3+4,UseFMA->False]=={Interval[{10,10}],Interval[{(6(-2 2^-53+2^-106)-4 2^-53)/10,(6(2 2^-53 +2^-106)+4 2^-53)/10}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[1+2]=={Interval[{3,3}],Interval[{-2^-53,2^-53}]}];
Assert[IEEEEvaluateWithRelativeError[2*3]=={Interval[{6,6}],Interval[{-2^-53,2^-53}]}];
Assert[IEEEEvaluateWithRelativeError[2/3]=={Interval[{CorrectlyRound[2/3],CorrectlyRound[2/3]}],Interval[{-2^-53,2^-53}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[-Interval[{-2,1}]]=={Interval[{-1,2}],Interval[{0,0}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[Interval[{-1,2}]^2]=={Interval[{0,4}],Interval[{-2^-53,2^-53}]}];
Assert[IEEEEvaluateWithRelativeError[Interval[{-1,3}]^3]=={Interval[{-1,27}],Interval[{-2 2^-53+2^-106,2 2^-53+2^-106}]}];
Assert[IEEEEvaluateWithRelativeError[Interval[{-1,3}]^4]=={Interval[{0,81}],Interval[{-2 2^-53+ 2^-106,2 2^-53+ 2^-106}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[Interval[{0,1}]+Interval[{0,2}]]==
	{Interval[{0,3}],Interval[{-2^-53,2^-53}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[Interval[{0,1}]+Interval[{-2,0}]]==
	{Interval[{-2,1}],Interval[{-\[Infinity],\[Infinity]}]}];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[
	{Interval[{5,5}],Interval[{-2^-53,2^-53}]}+
	{Interval[{6,6}],Interval[{-2^-53,2^-53}]}]==
	{Interval[{11,11}],Interval[{-2 2^-53+2^-106,2 2^-53+2^-106}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[CorrectlyRound[\[Pi]]]=={Interval[{884279719003555/281474976710656,884279719003555/281474976710656}],Interval[{0,0}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[
	Exact[
		{Interval[{5,5}],Interval[{-2^-53,2^-53}]}+
		{Interval[{6,6}],Interval[{-2^-53,2^-53}]}]]==
	{Interval[{11,11}],Interval[{-2^-53,2^-53}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[
	{Interval[{5,5}],Interval[{-2^-53,2^-53}]} *
	{Interval[{6,6}],Interval[{-2^-53,2^-53}]}]==
	{Interval[{30,30}],Interval[{-3 2^-53+3 2^-106-2^-159,3 2^-53+3 2^-106+2^-159}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[
Exact[
		{Interval[{5,5}],Interval[{-2^-53,2^-53}]}*
		{Interval[{6,6}],Interval[{-2^-53,2^-53}]}]]==
	{Interval[{30,30}],Interval[{-2 2^-53+2^-106,2 2^-53+2^-106}]}]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[\[Pi]]==$Failed];
Assert[IEEEEvaluateWithRelativeError[33333333333333333]==$Failed];


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[0.1`20]=={Interval[{3602879701896397/36028797018963968,3602879701896397/36028797018963968}],Interval[{(3602879701896397/36028797018963968)/0.1`20 -1,(3602879701896397/36028797018963968)/0.1`20 -1}]}]


(* ::Input::Initialization:: *)
Module[
{x=3},
Assert[IEEEEvaluateWithRelativeError[x]=={Interval[{3,3}],Interval[{0,0}]}]];


(* ::Input::Initialization:: *)
Module[
{a,b,poly,p,q,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly[u_]:=u^3(a+b u^2);
p=CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
	CorrectlyRound[
	CorrectlyRound[a]+
	CorrectlyRound[b]CorrectlyRound[w^2]]];
q=CorrectlyRound[poly[v]];
Assert[IEEEEvaluateWithRelativeError[poly[v]][[1]]!=Interval[{q,q}]];
Assert[IEEEEvaluateWithRelativeError[poly[v]][[1]]==Interval[{p,p}]]]


(* ::Input::Initialization:: *)
Module[
{a,b,poly,p,q,v,w},
a=-0.166666666666666666666648443029;
b=0.00833333316414932033147740536488;
v=0.000975`30;
w=CorrectlyRound[v];
poly=Function[u,u^3(a+b u^2)];
p=CorrectlyRound[CorrectlyRound[w CorrectlyRound[w^2]]
	CorrectlyRound[
	CorrectlyRound[a]+
	CorrectlyRound[b]CorrectlyRound[w^2]]];
q=CorrectlyRound[poly[v]];
Assert[IEEEEvaluateWithRelativeError[poly[v]][[1]]!=Interval[{q,q}]];
Assert[IEEEEvaluateWithRelativeError[poly[v]][[1]]==Interval[{p,p}]]]


(* ::Input::Initialization:: *)
Assert[IEEEEvaluateWithRelativeError[2,3]==$Failed];


(* ::Subsection:: *)
(*Higham' s \[Gamma] Model*)


(* ::Text:: *)
(*[Hig02], Lemma 3.1.*)


(* ::Input::Initialization:: *)
\[Gamma][n_]:=n 2^-53/(1-n 2^-53)


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}]Interval[{RandomReal[],RandomReal[]}]]][[2]]]<=\[Gamma][1],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[(Interval[{RandomReal[],RandomReal[]}]Interval[{RandomReal[],RandomReal[]}]) Interval[{RandomReal[],RandomReal[]}]]][[2]]]<=\[Gamma][2],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[((Interval[{RandomReal[],RandomReal[]}]Interval[{RandomReal[],RandomReal[]}]) Interval[{RandomReal[],RandomReal[]}])Interval[{RandomReal[],RandomReal[]}]]][[2]]]<=\[Gamma][3],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}]^2]][[2]]]<=\[Gamma][1],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}]^3]][[2]]]<=\[Gamma][2],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}]^4]][[2]]]<=\[Gamma][3],{100}],TrueQ]]


(* ::Text:: *)
(*[Hig02], equation 3.5.*)


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[
Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]][[2]]]]<=\[Gamma][2],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[
Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+(Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}])][[2]]]]<=\[Gamma][3],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[
Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+(Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}])][[2]]]]<=\[Gamma][3],{100}],TrueQ]]


(* ::Input::Initialization:: *)
Assert[AllTrue[Table[
Max[Abs[IEEEEvaluateWithRelativeError[Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+(Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+(Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]+Interval[{RandomReal[],RandomReal[]}] Interval[{RandomReal[],RandomReal[]}]))][[2]]]]<=\[Gamma][4],{100}],TrueQ]]
