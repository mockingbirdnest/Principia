(* ::Package:: *)

(* ::Title:: *)
(*Accurate Table Printing*)


(* ::Input:: *)
(*SetDirectory[ParentDirectory[NotebookDirectory[]]]*)


(* ::Input:: *)
(*Begin["sincos`"]*)


(* ::Input:: *)
(*<< "functions\\sin_cos_18.wl"*)


(* ::Input:: *)
(*End[]*)


(* ::Text:: *)
(*Returns a string representation of the 100-digit approximation of the given rational, using the C++ syntax for hex literals.*)


(* ::Input:: *)
(*hexFloatLiteral[x_Rational,signed_]:=*)
(* With[*)
(*  {groups=5,*)
(*group=5,*)
(*i=IntegerPart[N[x,100]],*)
(*f=FractionalPart[N[x,100]]},*)
(*  StringJoin[If[signed,If[x<0,"-","+"],""],"0x",IntegerString[i,16],".",StringRiffle[StringPartition[IntegerString[Round[f*16^(group*groups)],16,group*groups],UpTo[groups]],"'"],"p0"]]*)


(* ::Text:: *)
(*Prints the C++ header defining the accurate tables.*)


(* ::Input:: *)
(*Export["numerics\\accurate_tables.mathematica.h", *)
(*Module[{indices=Flatten[Position[ListQ/@Table[sincos`accurateTables[i],{i,0,500}],True]]-1,min,max,width},min=Min[indices];max=Max[indices];width=Ceiling[Log10[max]];*)
(*"#pragma once*)
(**)
(*#include <array>*)
(*#include <limits>*)
(**)
(*namespace principia {*)
(*namespace numerics {*)
(*namespace _accurate_tables {*)
(*namespace internal {*)
(**)
(*struct SinCosAccurateValues {*)
(*  double x;*)
(*  double sin_x;*)
(*  double cos_x;*)
(*};*)
(**)
(*constexpr std::array<SinCosAccurateValues, " <> ToString[max + 1] <> "> SinCosAccurateTable{{\n" <>*)
(*StringRiffle[Join[Table["    /*"<>StringPadLeft[ToString[i],width]<>"*/{.x =     std::numeric_limits<double>::signaling_NaN(),\n"<>StringRepeat[" ",width+9]<>".sin_x = std::numeric_limits<double>::signaling_NaN(),\n"<>StringRepeat[" ",width+9]<>".cos_x = std::numeric_limits<double>::signaling_NaN()",{i,0,min-1}],Table["    /*"<>StringPadLeft[ToString[i],width]<>"*/{.x =     " <> hexFloatLiteral[sincos`accurateTables[i][[1]],False] <> ",\n"<>StringRepeat[" ",width+9]<>".sin_x = " <> hexFloatLiteral[sincos`accurateTables[i][[2]],False] <> ",\n"<>StringRepeat[" ",width+9]<>".cos_x = "<> hexFloatLiteral[sincos`accurateTables[i][[3]],False], {i,indices}]],"},\n"]<>"}}};*)
(**)
(*}  // namespace internal*)
(**)
(*using internal::SinCosAccurateTable;*)
(**)
(*}  // namespace _accurate_tables*)
(*}  // namespace numerics*)
(*}  // namespace principia*)
(*"], "text"]*)
