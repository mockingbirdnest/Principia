(* ::Package:: *)

(* ::Section:: *)
(*Fukushima "Precise and fast computation of the general complete elliptic integral of the second kind", section 2.2*)


(* ::Input:: *)
(*fukushimaB[\[CurlyPhi]_,m_]:=(EllipticE[\[CurlyPhi],m]-(1-m) EllipticF[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaD[\[CurlyPhi]_,m_]:=(EllipticF[\[CurlyPhi],m]-EllipticE[\[CurlyPhi],m])/m*)


(* ::Input:: *)
(*fukushimaB[m_]:=fukushimaB[\[Pi]/2,m]*)


(* ::Input:: *)
(*fukushimaD[m_]:=fukushimaD[\[Pi]/2,m]*)


(* ::Input:: *)
(*fukushimaBstar[m_]:=m fukushimaB[m]*)


(* ::Input:: *)
(*fukushimaDstar[m_]:=m fukushimaD[m]*)


(* ::Input:: *)
(*fukushimaK0[mc_]:=(EllipticK[mc]/\[Pi])(-Log[16 EllipticNomeQ[mc]/mc])*)


(* ::Input:: *)
(*fukushimaKX[mc_]:=EllipticK[mc]/\[Pi]*)


(* ::Input:: *)
(*fukushimaX[mc_]:=-Log[mc/16]*)


(* ::Input:: *)
(*fukushimaE0[mc_]:=(\[Pi]/(2 EllipticK[mc]))+(1-(EllipticE[mc]/EllipticK[mc]))fukushimaK0[mc]*)


(* ::Input:: *)
(*fukushimaEX[mc_]:=(1-(EllipticE[mc]/EllipticK[mc]))fukushimaKX[mc]*)


(* ::Input:: *)
(*fukushimaBstar0[mc_]:=fukushimaE0[mc]-mc fukushimaK0[mc]*)


(* ::Input:: *)
(*fukushimaBstarX[mc_]:=fukushimaEX[mc]-mc fukushimaKX[mc]*)


(* ::Input:: *)
(*fukushimaDstar0[mc_]:=fukushimaK0[mc]-fukushimaE0[mc]*)


(* ::Input:: *)
(*fukushimaDstarX[mc_]:=fukushimaKX[mc]-fukushimaEX[mc]*)


(* ::Input:: *)
(*Series[fukushimaDstar0[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[fukushimaBstar0[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[fukushimaBstarX[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[fukushimaK0[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[fukushimaKX[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[fukushimaEX[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[EllipticNomeQ[mc],{mc,0,5}]*)


(* ::Input:: *)
(*Series[(EllipticK[mc]-1)/(\[Pi]/2),{mc,0,5}]*)


(* ::Input:: *)
(*Series[(EllipticK[mc]-EllipticE[mc])/(\[Pi]/2),{mc,0,5}]*)


(* ::Input:: *)
(*fukushimaEX[mc]//Simplify*)


(* ::Input:: *)
(*fukushimaKX[mc] // Simplify*)


(* ::Input:: *)
(*fukushimaBstarX[mc]//Simplify*)


(* ::Input:: *)
(*CoefficientList[Series[fukushimaKX[mc],{mc,0,7}],mc]*)


(* ::Input:: *)
(*CoefficientList[Series[fukushimaEX[mc],{mc,0,8}],mc]*)


(* ::Input:: *)
(*CoefficientList[Series[fukushimaBstarX[mc],{mc,0,8}],mc]*)
