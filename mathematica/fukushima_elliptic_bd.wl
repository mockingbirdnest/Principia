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
(*fukushimaBm1[\[Phi]_]:=Limit[Series[fukushimaB[\[Phi],mm],{mm,1,4},Assumptions->{mm\[Element]Reals,mm<1}],mm->1]*)


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


(* ::Section:: *)
(*Fukushima "Precise and fast computation of a general incomplete elliptic integral of the third kind by half and double argument transformations"*)


(* ::Input:: *)
(*bulirschCel[kc_,p_,a_,b_]:=Integrate[((a Cos[\[Theta]]^2+b Sin[\[Theta]]^2)/(Cos[\[Theta]]^2+p Sin[\[Theta]]^2))/Sqrt[Cos[\[Theta]]^2+kc^2 Sin[\[Theta]]^2],{\[Theta],0,\[Pi]/2}]*)


(* ::Input:: *)
(*fukushimaJ[\[Phi]_,n_,m_]:=(EllipticPi[n,\[Phi],m]-EllipticF[\[Phi],m])/n*)


(* ::Input:: *)
(*Series[bulirschCel[kc,nc,0,1],{kc,0,1},Assumptions->{kc\[Element]Reals,kc>0,nc\[Element]Reals,nc>0,nc<=1}]*)
