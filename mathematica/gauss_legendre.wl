(* ::Package:: *)

(* ::Text:: *)
(*See https://dlmf.nist.gov/3.5#v.*)


ClearAll[LegendreNodes];
LegendreNodes[n_]:=LegendreNodes[n]=Sort[N[Last/@Flatten[Solve[LegendreP[n,x]==0,x]],300]]


MonicLegendreP[n_][x_]:=LegendreP[n,x]/Coefficient[LegendreP[n,x],x^n]


ClearAll[LegendreWeights];
LegendreWeights[n_]:=LegendreWeights[n]=Table[NIntegrate[MonicLegendreP[n][x]/((x-LegendreNodes[n][[k]])MonicLegendreP[n]'[LegendreNodes[n][[k]]]),{x,-1,1},AccuracyGoal->\[Infinity],PrecisionGoal->100,WorkingPrecision->200,Method->{"GlobalAdaptive",Method->"GaussKronrodRule"}],{k,1,n}]


nodes=Table[
Table[
{n,k,Re[N[LegendreNodes[n][[k+1]],46]]},
{k,0,n-1}],
{n,0,50}];


nodes[[;;,;;,-1]]//TableForm


And@@Flatten[Table[nodes[[n+1,k+1,3]]==-nodes[[n+1,-(k+1),3]],{n,0,50},{k,0,n-1}]]


weights=Table[
Table[
{n,k,Re[N[LegendreWeights[n][[k+1]],46]]},
{k,0,n-1}],
{n,0,50}];


weights[[;;,;;,-1]]//TableForm


And@@Flatten[Table[weights[[n+1,k+1,3]]==weights[[n+1,-(k+1),3]],{n,0,50},{k,0,n-1}]]


decimalFloatLiteral[x_Real|x:0,exponentWidth_Integer,signed_]:=
 With[
  {m=MantissaExponent[x][[1]]*10,
   e=If[x==0,0,MantissaExponent[x][[2]]-1]},
  StringJoin[
   StringRiffle[#,"'"]&/@
    {{If[signed&&m>=0,"+",Nothing]},{#[[1]]},
     If[Length[#]>1,{"."},Nothing],
     If[Length[#]>1,StringPartition[#[[2]],UpTo[5]],Nothing],
     If[e!=0,{"e"<>If[e>0,"+","-"]<>IntegerString[e,10,exponentWidth]},Nothing]}&[
     StringSplit[ToString[m],"."]]]]


SetDirectory[NotebookDirectory[]]


Export[
"..\\numerics\\legendre_roots.mathematica.h",
"
#pragma once

#include \"numerics/fixed_arrays.hpp\"

namespace principia {
namespace numerics {
namespace _legendre_roots {
namespace internal {

using namespace principia::numerics::_fixed_arrays;

// Roots of the Legendre polynomials of degree up to 50.
constexpr FixedStrictlyLowerTriangularMatrix<double, "<>ToString[51]<>">
LegendreRoots{{{
"<>Flatten@Map[
With[
 {n=#[[1]],k=#[[2]],z=#[[3]]},
 "    /*"<>If[k==0,"n="<>StringPadLeft[ToString[n],2]<>", ","      "]<>"k="<>StringPadLeft[ToString[k],2]<>"*/"<>decimalFloatLiteral[z,1,True]<>",\n"]&,
nodes,{2}]<>"}}};

}  // namespace internal

using internal::LegendreRoots;

}  // namespace _legendre_roots
}  // namespace numerics
}  // namespace principia
",
"text"]


Export[
"..\\numerics\\gauss_legendre_weights.mathematica.h",
"
#pragma once

#include \"numerics/fixed_arrays.hpp\"

namespace principia {
namespace numerics {
namespace _gauss_legendre_weights {
namespace internal {

using namespace principia::numerics::_fixed_arrays;

// Weights for Gauss-Legendre quadrature with up to 50 points.
constexpr FixedStrictlyLowerTriangularMatrix<double, "<>ToString[51]<>">
GaussLegendreWeights{{{
"<>Flatten@Map[
With[
 {n=#[[1]],k=#[[2]],z=#[[3]]},
 "    /*"<>If[k==0,"n="<>StringPadLeft[ToString[n],2]<>", ","      "]<>"k="<>StringPadLeft[ToString[k],2]<>"*/"<>decimalFloatLiteral[z,1,False]<>",\n"]&,
weights,{2}]<>"}}};

}  // namespace internal

using internal::GaussLegendreWeights;

}  // namespace _gauss_legendre_weights
}  // namespace numerics
}  // namespace principia
",
"text"]
