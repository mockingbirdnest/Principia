(* ::Package:: *)

NormalizationFactor[n_,m_]:=Sqrt[((n-m)!(2n+1)(2-KroneckerDelta[0,m]))/(n+m)!]


pnrm[n_,m_,sin\[Phi]_]:=NormalizationFactor[n,m]LegendreP[n,m,sin\[Phi]]


N[Table[Table[Maximize[{Abs[pnrm[n,m,z]],z>=-1,z<=1},z][[1]],{m,0,n}],{n,0,5}],51]//TableForm


maxPnrm=Table[
Table[
{n,m,
SetPrecision[
Check[
NMaximize[
 {Abs[pnrm[n,m,z]],z>=-1,z<=1},
 {z,-1,1},
 PrecisionGoal->51,
 AccuracyGoal->\[Infinity],
 (* For WorkingPrecision\[Rule]103, all values converge
 except (39,1).  For WorkingPrecision\[Rule]104, many
 values fail to converge.  Don't ask. *)
 WorkingPrecision->If[n==39&&m==1,104,103]][[1]],
ToString[{n,m}]<>": error"],
51]},
{m,0,n}],
{n,0,50}];
Map[Last,maxPnrm,{2}]//TableForm


decimalFloatLiteral[x_Real]:=
 With[
  {m=MantissaExponent[x][[1]]*10,
   e=MantissaExponent[x][[2]]-1},
  StringJoin[
   StringRiffle[#,"'"]&/@
    {{#[[1]]},
     If[Length[#]>1,{"."},Nothing],
     If[Length[#]>1,StringPartition[#[[2]],UpTo[5]],Nothing],
     If[e!=0,"e"<>ToString[e],Nothing]}&[
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
 "    /*"<>If[m==0,"n="<>StringPadLeft[ToString[n],2]<>", ","      "]<>"m="<>StringPadLeft[ToString[m],2]<>"*/"<>decimalFloatLiteral[z]<>",\n"]&,
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
       decimalFloatLiteral[N[NormalizationFactor[n,m],46]]<>",\n",
   {m,0,n}],
  {n,0,50}]]<>"}}};

}  // namespace numerics
}  // namespace principia
",
"text"]
