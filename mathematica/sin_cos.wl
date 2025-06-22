(* ::Package:: *)

(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]];Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point_evaluation.wl"}]];*)


(* ::Input:: *)
(*SetFloatingPointFormat[binary64];*)
(*SetRoundingMode[NearestTiesToEven];*)


(* ::Input:: *)
(*On[Assert]*)


(* ::Input:: *)
(*?IEEE754FloatingPoint`**)


(* ::Input:: *)
(*?UlpDistance*)


(* ::Input:: *)
(*M=binary64[[1]]*)


(* ::Input:: *)
(*\[GothicU][x_]:=FromRepresentation[Representation[x]+1]-x*)


(* ::Input:: *)
(*Assert[\[GothicU][1]==2^(1-M)]*)


(* ::Section:: *)
(*Approximation of \[Pi]/2*)


(* ::Subsection:: *)
(*Two - Term Approximation*)


(* ::Input:: *)
(*Bits[\[Pi]/2]*)


(* ::Input:: *)
(*\[Kappa]1=8;*)
(*\[Kappa]\[Prime]1=5;*)


(* ::Input:: *)
(*C1=Truncate[\[Kappa]1,\[Pi]/2];*)


(* ::Input:: *)
(*Assert[C1==Truncate[\[Kappa]\[Prime]1,\[Pi]/2]]*)


(* ::Input:: *)
(*\[Delta]C1=CorrectlyRound[\[Pi]/2-C1];*)


(* ::Input:: *)
(*Assert[Abs[\[Pi]/2-C1-\[Delta]C1]<=2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]]*)


(* ::Text:: *)
(*This is very pessimistic:*)


(* ::Input:: *)
(*N[{Abs[\[Pi]/2-C1-\[Delta]C1],2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C1]<2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C1],2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*HexLiteral[C1]*)


(* ::Input:: *)
(*HexLiteral[\[Delta]C1]*)


(* ::Subsection:: *)
(*Three - Term Approximation*)


(* ::Input:: *)
(*Bits[\[Pi]/2]*)


(* ::Text:: *)
(*The paper has wrong values for \[Kappa]\[Prime]2 and \[Kappa]\[DoublePrime]2:*)


(* ::Input:: *)
(*\[Kappa]2=18;*)
(*\[Kappa]\[Prime]2=14;*)
(*\[Kappa]\[DoublePrime]2=15;*)


(* ::Input:: *)
(*C2=Truncate[\[Kappa]2,\[Pi]/2];*)


(* ::Input:: *)
(*Assert[C2==Truncate[\[Kappa]\[Prime]2,\[Pi]/2]]*)


(* ::Input:: *)
(*C\[Prime]2=Truncate[\[Kappa]2,\[Pi]/2-C2];*)


(* ::Input:: *)
(*Assert[C\[Prime]2==Truncate[\[Kappa]\[DoublePrime]2,\[Pi]/2-C2]]*)


(* ::Input:: *)
(*\[Delta]C2=CorrectlyRound[\[Pi]/2-C2-C\[Prime]2];*)


(* ::Input:: *)
(*Assert[Abs[\[Pi]/2-C2-C\[Prime]2-\[Delta]C2]<=2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2 M-1)\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Pi]/2-C2-C\[Prime]2-\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2 M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C2]<2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-53)(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-53)(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*HexLiteral[C2]*)


(* ::Input:: *)
(*HexLiteral[C\[Prime]2]*)


(* ::Input:: *)
(*HexLiteral[\[Delta]C2]*)
