(* ::Package:: *)

NormalizationFactor[n_,m_]:=Sqrt[((n-m)!(2n+1)(2-KroneckerDelta[0,m]))/(n+m)!]


pnrm[n_,m_,sin\[Phi]_]:=NormalizationFactor[n,m]LegendreP[n,m,sin\[Phi]]


(maximizations=Table[Table[Maximize[{Abs[pnrm[n,m,z]],z>=-1,z<=1},z][[1]],{m,0,n}],{n,0,5}]);
N[maximizations,51]//TableForm


Show[
Plot[Evaluate@Table[pnrm[#,m,z],{m,0,#}],{z,-1,1}],
Plot[Evaluate[-maximizations[[#+1,;;]]],{z,-1,1}],ImageSize->500]&/@{4,5}//Row


ClearAll[maxP];
maxP[n_,n_]:=maxP[n,n]={Abs[pnrm[n,n,0]],0};
maxP[n_,m_]:=maxP[n,m]=Check[
{#[[1]],#[[2,1,2]]}&@NMaximize[
 {Abs[pnrm[n,m,z]],z>=-1,z<=SetPrecision[maxP[n,m+1][[2]],\[Infinity]]},
 {z,-1,maxP[n,m+1][[2]]},
 PrecisionGoal->51,
 AccuracyGoal->\[Infinity],
 WorkingPrecision->104],
ToString[{n,m}]<>": error"]


(nmaximizations=Table[Table[maxP[n,m][[1]],{m,0,n}],{n,0,5}])//TableForm


maximizations-nmaximizations//TableForm


Table[
ToExpression[StringReplace[ToString[N[nmaximizations,sigdec]],"."->""]]/10^(sigdec-1)-Round[maximizations,10^-(sigdec-1)],
{sigdec,{46,51,103}}]//Column


maxPnrm=Table[
Table[
{n,m,N[maxP[n,m][[1]],46]},
{m,0,n}],
{n,0,50}];
Map[Last,maxPnrm,{2}]//TableForm


decimalFloatLiteral[x_Real,exponentWidth_Integer]:=
 With[
  {m=MantissaExponent[x][[1]]*10,
   e=MantissaExponent[x][[2]]-1},
  StringJoin[
   StringRiffle[#,"'"]&/@
    {{#[[1]]},
     If[Length[#]>1,{"."},Nothing],
     If[Length[#]>1,StringPartition[#[[2]],UpTo[5]],Nothing],
     If[e!=0,"e"<>If[e>0,"+","-"]<>IntegerString[e,10,exponentWidth],Nothing]}&[
     StringSplit[ToString[m],"."]]]]


SetDirectory[NotebookDirectory[]]


Export[
"..\\numerics\\max_abs_normalized_associated_legendre_function.mathematica.h",
"
#pragma once

#include \"numerics/fixed_arrays.hpp\"

namespace principia {
namespace numerics {

// Global maxima over [-1, 1] of the absolute value of the normalized associated
// Legendre functions.
constexpr FixedLowerTriangularMatrix<double, "<>ToString[51]<>">
MaxAbsNormalizedAssociatedLegendreFunction{{{
"<>Flatten@Map[
With[
 {n=#[[1]],m=#[[2]],z=#[[3]]},
 "    /*"<>If[m==0,"n="<>StringPadLeft[ToString[n],2]<>", ","      "]<>"m="<>StringPadLeft[ToString[m],2]<>"*/"<>decimalFloatLiteral[z,1]<>",\n"]&,
maxPnrm,{2}]<>"}}};

}  // namespace numerics
}  // namespace principia
",
"text"]


Export[
"..\\numerics\\legendre_normalization_factor.mathematica.h",
"
#pragma once

#include \"numerics/fixed_arrays.hpp\"

namespace principia {
namespace numerics {

// Multiplying a normalized Cnm or Snm coefficient by this factor yields an
// unnormalized coefficient.  Dividing an unnormalized Cnm or Snm coefficient by
// this factor yields a normalized coefficient.
constexpr FixedLowerTriangularMatrix<double, "<>ToString[51]<>">
LegendreNormalizationFactor{{{
"<>Flatten[
 Table[
  Table[
   "    /*"<>If[m==0,"n="<>StringPadLeft[ToString[n],2]<>", ","      "]<>"m="<>StringPadLeft[ToString[m],2]<>"*/"<>
       decimalFloatLiteral[N[NormalizationFactor[n,m],46],2]<>",\n",
   {m,0,n}],
  {n,0,50}]]<>"}}};

}  // namespace numerics
}  // namespace principia
",
"text"]
