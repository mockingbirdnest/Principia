(* ::Package:: *)

(* ::Input:: *)
(*Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]];*)


(* ::Input:: *)
(*SetFloatingPointFormat[binary64]*)


(* ::Input:: *)
(*SetRoundingMode[NearestTiesToEven]*)


(* ::Text:: *)
(*[Mul16] section 2.1.3*)


(* ::Input:: *)
(*mullerULP[x_,fmt_]:=Module[{e=Floor[Log2[x]]},2^(e-fmt[[1]]+1)]*)


(* ::Text:: *)
(*[Mul16] section 11.2.2*)


(* ::Input:: *)
(*bdlC = \[Pi]/2*)


(* ::Input:: *)
(*bdlR=CorrectlyRound[1/bdlC]*)


(* ::Input:: *)
(*HexLiteral[bdlR,Quotes->4]*)


(* ::Input:: *)
(*bdlC1=Round[1/(4 bdlR  mullerULP[1/bdlR,binary64])]4 mullerULP[1/bdlR,binary64]*)


(* ::Input:: *)
(*HexLiteral[bdlC1,Quotes->4]*)


(* ::Input:: *)
(*bdlC2=Round[(bdlC-bdlC1)/(2^(-binary64[[1]]+4)mullerULP[bdlC1,binary64])]2^(-binary64[[1]]+4)mullerULP[bdlC1,binary64]*)


(* ::Input:: *)
(*HexLiteral[bdlC2,Quotes->4]*)


(* ::Input:: *)
(*bdlC3=CorrectlyRound[bdlC-(bdlC1+bdlC2)]*)


(* ::Input:: *)
(*HexLiteral[bdlC3,Quotes->4]*)


(* ::Input:: *)
(*Log2[mullerULP[bdlC3,binary64]]//N*)


(* ::Input:: *)
(*Log2[bdlC3]//N*)


(* ::Input:: *)
(*bdlAddend=3 2^(binary64[[1]]-2)*)


(* ::Input:: *)
(*HexLiteral[bdlAddend,Quotes->4]*)


(* ::Input:: *)
(*bdl\[Epsilon]=2^(-binary64[[1]]-18)*)


(* ::Input:: *)
(*bdlThreshold=CorrectlyRound[(2^( binary64[[1]]-3) bdl\[Epsilon]/mullerULP[bdlC,binary64]-1)bdlC,RoundingMode->Toward0]*)


(* ::Input:: *)
(*HexLiteral[bdlThreshold,Quotes->4]*)


(* ::Input:: *)
(*\[Kappa]3=18*)


(* ::Input:: *)
(*bdl2Multiplier =2^(\[Kappa]3-binary64[[1]]+4)*)


(* ::Input:: *)
(*Log2[bdl2Multiplier]*)


(* ::Input:: *)
(*N[bdl2Multiplier]*)


(* ::Input:: *)
(*bdl3Multiplier=2^(\[Kappa]3-2 binary64[[1]]+4)(3 + 2^(\[Kappa]3+1))*)


(* ::Input:: *)
(*Log2[Denominator[bdl3Multiplier]]*)


(* ::Input:: *)
(*N[bdl3Multiplier]*)


(* ::Input:: *)
(*bdlThreshold=\[Pi]^2 2^(binary64[[1]]-\[Kappa]3-7)*)
