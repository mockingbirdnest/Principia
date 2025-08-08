(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


(* ::Input:: *)
(*<<FunctionApproximations`*)


(* ::Input:: *)
(*Get[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"functions","sin_cos_18_only1.wl"}]];Get[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"functions","sin_cos_18_not1.wl"}]];Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point.wl"}]];Get[FileNameJoin[{NotebookDirectory[],"ieee754_floating_point_evaluation.wl"}]];*)


(* ::Input:: *)
(*SetFloatingPointFormat[binary64];*)
(*SetRoundingMode[NearestTiesToEven];*)


(* ::Input:: *)
(*On[Assert]*)


(* ::Text:: *)
(*Extra bits for argument reduction:*)


(* ::Input:: *)
(*\[Kappa]3=18*)


(* ::Text:: *)
(*Characteristics of the floating-point model:*)


(* ::Input:: *)
(*M=binary64[[1]]*)


(* ::Input:: *)
(*\[GothicU][x_]:=FromRepresentation[Representation[x]+1]-x*)


(* ::Input:: *)
(*Assert[\[GothicU][1]==2^(1-M)]*)


(* ::Input:: *)
(*\[Gamma][n_]:=n 2^-M/(1-n 2^-M)*)


(* ::Input:: *)
(*uInterval:=Interval[{-\[GothicU][1/2],\[GothicU][1/2]}]*)


(* ::Input:: *)
(*binaryBounds[bounds_List]:=Module[*)
(*{b=N[bounds,20],min,max},*)
(*min=Min[b];max=Max[b];{If[min<0,"-",""]<>"2^"<> ToString[Log2[Abs[min]]],If[max<0,"-",""]<>"2^"<>ToString[Log2[Abs[max]]]}];*)
(*binaryBounds[bounds_Interval]:=binaryBounds[MinMax[bounds]];*)


(* ::Subsection::Closed:: *)
(*Accurate Tables Utilities*)


(* ::Input:: *)
(*hasAccurateTables[i_]:=ValueQ[accurateTables[i],Method->"TrialEvaluation"]*)


(* ::Input:: *)
(*accurateTablesMaxIndex=Do[If[!hasAccurateTables[i],Return[i-1]],{i,1,500}]*)


(* ::Input:: *)
(*accurateTablesStep=1/512;*)


(* ::Text:: *)
(*Find the bounds over which we'll have to evaluate the polynomial, based on the excursion away from multiples of 1/512:*)


(* ::Input:: *)
(*accurateTablesXIntervals=Table[Interval[{(2i-1)accurateTablesStep/2,(2i+1)accurateTablesStep/2}],{i,1,accurateTablesMaxIndex}];*)


(* ::Input:: *)
(*accurateTablesHIntervals=Table[accurateTablesXIntervals[[i]]-accurateTables[i][[1]],{i,1,accurateTablesMaxIndex}];*)


(* ::Input:: *)
(*{hLB,hUB}={Min[Min/@accurateTablesHIntervals],Max[Max/@accurateTablesHIntervals]};*)


(* ::Input:: *)
(*hInterval=Interval[{hLB,hUB}];*)


(* ::Input:: *)
(*N[Log2[Abs[{hLB,hUB}]]]*)


(* ::Input:: *)
(* N[Log2[Abs[{hLB,hUB}]-accurateTablesStep/2]]*)


(* ::Text:: *)
(*Check the Sterbenz condition for computing s0+c0 h exactly (the subtraction h' - s0 is exact):*)


(* ::Input:: *)
(*AllTrue[Table[Module[{*)
(*t=accurateTables[i],s0,c0,h,h\[Prime]},*)
(*s0=t[[2]];c0=t[[3]];h=accurateTablesHIntervals[[i]];h\[Prime]=(s0+c0 h)(1+uInterval);s0/2<=h\[Prime]<=2s0],{i,1,accurateTablesMaxIndex}],TrueQ]*)


(* ::Text:: *)
(*Similarly the Sterbenz condition for computing c0 - h s0 exactly (the subtraction h'- c0 is exact):*)


(* ::Input:: *)
(*AllTrue[Table[Module[{*)
(*t=accurateTables[i],s0,c0,h,h\[Prime]},*)
(*s0=t[[2]];*)
(*c0=t[[3]];*)
(*h=accurateTablesHIntervals[[i]];*)
(*h\[Prime]=(c0-h s0)(1+uInterval);*)
(*c0/2<=h\[Prime]<=2 c0],{i,1,accurateTablesMaxIndex}],TrueQ]*)


(* ::Section:: *)
(*Approximation of \[Pi]/2*)


(* ::Text:: *)
(*This section verifies the error bounds for the approximations of \[Pi]/2.*)


(* ::Subsection::Closed:: *)
(*Two - Term Approximation*)


(* ::Input:: *)
(*\[Kappa]1;\[Kappa]\[Prime]1;C1;\[Delta]C1;*)


(* ::Input:: *)
(*Begin["TwoTermApproximation`"]*)


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
(*The upper bound on the error of the approximation is very pessimistic:*)


(* ::Input:: *)
(*N[{Abs[\[Pi]/2-C1-\[Delta]C1],2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Text:: *)
(*That's probably because there are three zeroes after \[Delta]C1:*)


(* ::Input:: *)
(*Bits[\[Pi]/2-C1]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C1]<2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C1],2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*HexLiteral[C1,Quotes->4]*)


(* ::Input:: *)
(*HexLiteral[\[Delta]C1,Quotes->4]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsection::Closed:: *)
(*Three - Term Approximation*)


(* ::Input:: *)
(*\[Kappa]2;\[Kappa]\[Prime]2;\[Kappa]\[DoublePrime]2;C2;C\[Prime]2;\[Delta]C2;*)


(* ::Input:: *)
(*Begin["ThreeTermApproximation`"]*)


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


(* ::Text:: *)
(*Here the upper bound on the error of the approximation is tight:*)


(* ::Input:: *)
(*N[{Abs[\[Pi]/2-C2-C\[Prime]2-\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2 M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Text:: *)
(*That's because the first bit after \[Delta]C2 is one:*)


(* ::Input:: *)
(*Bits[\[Pi]/2-C2-C\[Prime]2]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C2]<2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M)(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M)(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*HexLiteral[C2,Quotes->4]*)


(* ::Input:: *)
(*HexLiteral[C\[Prime]2,Quotes->4]*)


(* ::Input:: *)
(*HexLiteral[\[Delta]C2,Quotes->4]*)


(* ::Input:: *)
(*End[]*)


(* ::Section:: *)
(*Argument Reduction*)


(* ::Text:: *)
(*This section verifies the error bounds on argument reduction.*)


(* ::Subsection::Closed:: *)
(*Argument Reduction Using the Two-Term Approximation*)


(* ::Input:: *)
(*Begin["TwoTermReduction`"]*)


(* ::Input:: *)
(*xInterval=Interval[{-2^\[Kappa]1 CorrectlyRound[\[Pi]/2],2^\[Kappa]1 CorrectlyRound[\[Pi]/2]}];*)


(* ::Text:: *)
(*Bound on n:*)


(* ::Input:: *)
(*nIntervals=IEEEEvaluateWithAbsoluteError[xInterval CorrectlyRound[2/\[Pi]]]*)


(* ::Input:: *)
(*nIntervalBeforeRounding=nIntervals[[1]]+nIntervals[[2]];*)


(* ::Input:: *)
(*Assert[Max[Abs[nIntervalBeforeRounding]]<=2^\[Kappa]1(1+\[Gamma][3])]*)


(* ::Input:: *)
(*N[{Max[Abs[nIntervalBeforeRounding]],2^\[Kappa]1(1+\[Gamma][3])},20]*)


(* ::Input:: *)
(*nInterval=Round[nIntervalBeforeRounding]*)


(* ::Text:: *)
(*The reason why the error bound on n is overly broad is that the rounding errors on the constants largely cancel:*)


(* ::Input:: *)
(*\[Delta]1=CorrectlyRound[\[Pi]/2]/(\[Pi]/2)-1;*)


(* ::Input:: *)
(*\[Delta]2=CorrectlyRound[2/\[Pi]]/(2/\[Pi])-1;*)


(* ::Input:: *)
(*N[{\[Delta]1,\[Delta]2}/\[GothicU][1/2],20]*)


(* ::Text:: *)
(*Bounds on misrounding.  Note that the second kind always has a larger error:*)


(* ::Input:: *)
(*misroundingBound1=(\[Pi]/4)(\[Gamma][2]/(1+\[Gamma][2]))(2^(\[Kappa]1+1)-1);*)


(* ::Input:: *)
(*N[misroundingBound1]*)


(* ::Input:: *)
(*misroundingBound2=(\[Pi]/4)\[Gamma][2](2^(\[Kappa]1+1)+1);*)


(* ::Input:: *)
(*N[misroundingBound2]*)


(* ::Text:: *)
(*Check that misrounding stays in the last interval:*)


(* ::Input:: *)
(*Assert[IntervalMemberQ[ accurateTablesXIntervals[[accurateTablesMaxIndex]],\[Pi]/4+misroundingBound1]]*)


(* ::Input:: *)
(*N[{ accurateTablesXIntervals[[accurateTablesMaxIndex]],\[Pi]/4+misroundingBound1}]*)


(* ::Text:: *)
(*Bounds on \[Delta]y:*)


(* ::Input:: *)
(*\[Delta]yInterval=IEEEEvaluateWithAbsoluteError[nInterval \[Delta]C1];*)


(* ::Input:: *)
(*Assert[Max[\[Delta]yInterval[[2]]]<=Max[\[Delta]yInterval[[1]]]\[Gamma][1]]*)


(* ::Input:: *)
(*N[{Max[\[Delta]yInterval[[2]]],Max[\[Delta]yInterval[[1]]]\[Gamma][1]}]*)


(* ::Text:: *)
(*Error on the approximation of \[Pi]/2:*)


(* ::Input:: *)
(*\[Zeta]=C1+\[Delta]C1-\[Pi]/2*)


(* ::Text:: *)
(*The bound on \[Zeta] is pessimistic because C1 + \[Delta]C1 is an unusually good approximation:*)


(* ::Input:: *)
(*Assert[\[Zeta]<=2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{\[Zeta],2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Text:: *)
(*Error on the overall reduction.  As expected, the bound is pessimistic:*)


(* ::Input:: *)
(*reductionInterval=\[Zeta] nInterval + \[Delta]yInterval[[2]];*)


(* ::Input:: *)
(*Assert[Max[reductionInterval]<2^(\[Kappa]1+\[Kappa]\[Prime]1-M+1)\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Max[reductionInterval],2^(\[Kappa]1+\[Kappa]\[Prime]1-M+1)\[GothicU][\[Pi]/2]},20]*)


(* ::Text:: *)
(*Threshold for dangerous rounding:*)


(* ::Input:: *)
(*angleReducedThreshold = 2^(\[Kappa]1+\[Kappa]\[Prime]1+\[Kappa]3-M+2)*)


(* ::Input:: *)
(*Denominator[angleReducedThreshold]//Log2*)


(* ::Input:: *)
(*N[angleReducedThreshold]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsection::Closed:: *)
(*Argument Reduction Using the Three-Term Approximation*)


(* ::Input:: *)
(*Begin["ThreeTermReduction`"]*)


(* ::Input:: *)
(*xInterval=Interval[{-2^\[Kappa]2 CorrectlyRound[\[Pi]/2],2^\[Kappa]2 CorrectlyRound[\[Pi]/2]}];*)


(* ::Text:: *)
(*Bound on n:*)


(* ::Input:: *)
(*nIntervals=IEEEEvaluateWithAbsoluteError[xInterval CorrectlyRound[2/\[Pi]]]*)


(* ::Input:: *)
(*nIntervalBeforeRounding=nIntervals[[1]]+nIntervals[[2]];*)


(* ::Input:: *)
(*Assert[Max[Abs[nIntervalBeforeRounding]]<=2^\[Kappa]2(1+\[Gamma][3])]*)


(* ::Input:: *)
(*N[{Max[Abs[nIntervalBeforeRounding]],2^\[Kappa]2(1+\[Gamma][3])},20]*)


(* ::Input:: *)
(*nInterval=Round[nIntervalBeforeRounding]*)


(* ::Text:: *)
(*Bounds on misrounding.  Note that the second kind always has a larger error:*)


(* ::Input:: *)
(*misroundingBound1=(\[Pi]/4)(\[Gamma][2]/(1+\[Gamma][2]))(2^(\[Kappa]2+1)-1);*)


(* ::Input:: *)
(*N[misroundingBound1]*)


(* ::Input:: *)
(*misroundingBound2=(\[Pi]/4)\[Gamma][2](2^(\[Kappa]2+1)+1);*)


(* ::Input:: *)
(*N[misroundingBound2]*)


(* ::Input:: *)
(*misroundingBound=Max[misroundingBound1,misroundingBound2];*)


(* ::Text:: *)
(*Bounds on \[Delta]y:*)


(* ::Input:: *)
(*\[Delta]yInterval=IEEEEvaluateWithAbsoluteError[nInterval \[Delta]C2];*)


(* ::Input:: *)
(*Assert[Max[\[Delta]yInterval[[2]]]<=Max[\[Delta]yInterval[[1]]]\[Gamma][1]]*)


(* ::Input:: *)
(*N[{Max[\[Delta]yInterval[[2]]],Max[\[Delta]yInterval[[1]]]\[Gamma][1]}]*)


(* ::Text:: *)
(*Error on the approximation of \[Pi]/2:*)


(* ::Input:: *)
(*\[Zeta]1=C2+C\[Prime]2+\[Delta]C2-\[Pi]/2*)


(* ::Input:: *)
(*Assert[\[Zeta]1<=2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2M-1)\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{\[Zeta]1,2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2M-1)\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*\[Zeta]2=2^(2-2M);*)


(* ::Text:: *)
(*The paper should use misrounding in the following bound!*)


(* ::Text:: *)
(*Error on the overall reduction:*)


(* ::Input:: *)
(*reductionInterval=(\[Zeta]1 nInterval + \[Delta]yInterval[[2]])(1+\[Zeta]2)+Interval[{-\[Pi]/4-misroundingBound,\[Pi]/4+misroundingBound}]\[Zeta]2;*)


(* ::Input:: *)
(*Assert[Max[reductionInterval]<2^(\[Kappa]2-2M)(2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2+1)+3)\[GothicU][\[Pi]/2]+2^(-2M)\[Pi]]*)


(* ::Input:: *)
(*N[{Max[reductionInterval],2^(\[Kappa]2-2M)(2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2+1)+3)\[GothicU][\[Pi]/2]+2^(-2M)\[Pi]},20]*)


(* ::Text:: *)
(*Threshold for dangerous rounding:*)


(* ::Input:: *)
(*angleReducedThreshold = 2^(\[Kappa]3-M)(2^(\[Kappa]2+\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M+2)+4)*)


(* ::Input:: *)
(*Denominator[angleReducedThreshold]//Log2*)


(* ::Input:: *)
(*N[angleReducedThreshold]*)


(* ::Input:: *)
(*End[]*)


(* ::Section:: *)
(*Polynomials*)


(* ::Text:: *)
(*This section derives the minimax polynomials used to approximate Sin and Cos.*)


(* ::Subsection::Closed:: *)
(*Helper Functions*)


(* ::Text:: *)
(*For Sin, we want to compute minimax polynomials that are odd and where the term in t has coefficient 1 so that the computation produces two terms, t plus a correction.  Therefore, we approximate the function sinFn below.  Note that it is singular near 0 and must therefore be replaced by its limit there.*)


(* ::Input:: *)
(*Series[Sin[t],{t,0,5}]*)


(* ::Input:: *)
(*sinFn[t_]:=(Sin[t]-t)/t^3*)


(* ::Input:: *)
(*Limit[sinFn[t],t->0]*)


(* ::Input:: *)
(*Series[sinFn[t],{t,0,2}]*)


(* ::Text:: *)
(*For Cos, we want to compute a minimax polynomial that is even and where the constant term is 1.  Therefore, we approximate the function cosFn below.  Note that it is singular near 0 and must therefore be replaced by its limit there.*)


(* ::Input:: *)
(*Series[Cos[t],{t,0,5}]*)


(* ::Input:: *)
(*cosFn[t_]:=(Cos[t]-1)/t^2*)


(* ::Input:: *)
(*Limit[cosFn[t],t->0]*)


(* ::Input:: *)
(*Series[cosFn[t],{t,0,2}]*)


(* ::Subsection::Closed:: *)
(*Sin Near  Zero*)


(* ::Text:: *)
(*Near zero the Gal and Bachelis method is not usable, so we use a plain minimax polynomial for the entire function.  We use the sinFn function with an error function suitable for the relative error.*)


(* ::Input:: *)
(*sin0ApproximationResult;*)
(*sin0Polynomial;*)
(*x0Interval;*)


(* ::Input:: *)
(*Begin["PolynomialSinNearZero`"]*)


(* ::Input:: *)
(*x0Max=1/1024;*)


(* ::Input:: *)
(*x0Min=-x0Max;*)


(* ::Input:: *)
(*x0Interval=Interval[{x0Min,x0Max}]*)


(* ::Input:: *)
(*sin0ApproximationResult=GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1,Sin[t]/t^3]},*)
(*{t,{0,x0Max},1,0},*)
(*x,WorkingPrecision->30]*)


(* ::Input:: *)
(*sin0Polynomial=Function[u, Evaluate[sin0ApproximationResult[[2,1]]/.x->u]]*)


(* ::Input:: *)
(*sin0ExactRelativeError = Function[u,1-(u+ u^3 sin0Polynomial[u^2])/Sin[u]]*)


(* ::Input:: *)
(*sin0IEEERelativeError = Function[u,1-(u+IEEEEvaluate[ u^3 sin0Polynomial[u^2]])/Sin[u]]*)


(* ::Text:: *)
(*The plot of the relative error with floating-point evaluation is very noisy and there is no reliable way to find its extrema:*)


(* ::Input:: *)
(*Plot[{sin0ExactRelativeError[x],sin0IEEERelativeError[x]},{x,x0Min,x0Max},WorkingPrecision->30,PlotPoints->200,PlotRange->Full]*)


(* ::Input:: *)
(*sin0RoundedCoefficients=Map[CorrectlyRound,CoefficientList[sin0Polynomial[x],x]]*)


(* ::Input:: *)
(*s3=sin0RoundedCoefficients[[1]];*)


(* ::Input:: *)
(*HexLiteral[s3,Quotes->4]*)


(* ::Input:: *)
(*s5=sin0RoundedCoefficients[[2]];*)


(* ::Input:: *)
(*HexLiteral[s5,Quotes->4]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsection::Closed:: *)
(*Sin Around Table Entries*)


(* ::Text:: *)
(*We we use the sinFn function with an error function suitable for the absolute error.*)


(* ::Input:: *)
(*cosApproximationResult;*)
(*cosPolynomial;*)
(*sinApproximationResult;*)
(*sinPolynomial;*)


(* ::Input:: *)
(*Begin["PolynomialSinAroundTableEntries`"]*)


(* ::Text:: *)
(*Build a symmetrical interval, because our polynomial is going to be odd anyway:*)


(* ::Input:: *)
(*hMax=Max[Abs[{hLB,hUB}]];*)


(* ::Input:: *)
(*hMin=-hMax;*)


(* ::Input:: *)
(*sinApproximationResult=GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,-1/6,sinFn[t]]},*)
(*{t,{0,hMax},1,0},*)
(*x,WorkingPrecision->30]*)


(* ::Input:: *)
(*sinPolynomial=Function[u, Evaluate[ sinApproximationResult[[2,1]]/.x->u]]*)


(* ::Input:: *)
(*sinExactRelativeError = Function[u,(u+ u^3 sinPolynomial[u^2])/Sin[u]-1]*)


(* ::Input:: *)
(*sinIEEERelativeError = Function[u,(u+IEEEEvaluate[ u^3 sinPolynomial[u^2]])/Sin[u]-1]*)


(* ::Text:: *)
(*Despite the funky error function, the relative error on the overall Sin is reasonable:*)


(* ::Input:: *)
(*Plot[{sinExactRelativeError[x],sinIEEERelativeError[x]},{x,hMin,hMax},WorkingPrecision->30,PlotPoints->200,PlotRange->Full]*)


(* ::Input:: *)
(*sinRoundedCoefficients=Map[CorrectlyRound,CoefficientList[sinPolynomial[x],x]]*)


(* ::Input:: *)
(*s3=sinRoundedCoefficients[[1]];*)


(* ::Input:: *)
(*HexLiteral[s3,Quotes->4]*)


(* ::Input:: *)
(*s5=sinRoundedCoefficients[[2]];*)


(* ::Input:: *)
(*HexLiteral[s5,Quotes->4]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsection::Closed:: *)
(*Cos Around Table Entries*)


(* ::Text:: *)
(*We use the cosFn function with an error function suitable for the absolute error.*)


(* ::Input:: *)
(*Begin["PolynomialCosAroundTableEntries`"]*)


(* ::Text:: *)
(*Build a symmetrical interval, because our polynomial is going to be even anyway:*)


(* ::Input:: *)
(*hMax=Max[Abs[{hLB,hUB}]];*)


(* ::Input:: *)
(*hMin=-hMax;*)


(* ::Input:: *)
(*cosApproximationResult=GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/2,cosFn[t]],If[t==0,-1/2,cosFn[t]]},*)
(*{t,{0,hMax},1,0},*)
(*x,WorkingPrecision->30]*)


(* ::Input:: *)
(*cosPolynomial=Function[u, Evaluate[ cosApproximationResult[[2,1]]/.x->u]]*)


(* ::Input:: *)
(*cosExactRelativeError = Function[u,(1+ u^2cosPolynomial[u^2])/Cos[u]-1]*)


(* ::Input:: *)
(*cosIEEERelativeError = Function[u,(1+IEEEEvaluate[ u^2 cosPolynomial[u^2]])/Cos[u]-1]*)


(* ::Text:: *)
(*Despite the funky error function, the relative error on the overall Cos is reasonable:*)


(* ::Input:: *)
(*Plot[{cosExactRelativeError[x],cosIEEERelativeError[x]},{x,hMin,hMax},WorkingPrecision->30,PlotPoints->200,PlotRange->Full]*)


(* ::Input:: *)
(*cosRoundedCoefficients=Map[CorrectlyRound,CoefficientList[cosPolynomial[x],x]]*)


(* ::Text:: *)
(*Note that the fact that this coefficient is a power of two has no effect on the accuracy as it's used as part of an FMA.*)


(* ::Input:: *)
(*c2=cosRoundedCoefficients[[1]];*)


(* ::Input:: *)
(*HexLiteral[c2,Quotes->4]*)


(* ::Input:: *)
(*c4=cosRoundedCoefficients[[2]];*)


(* ::Input:: *)
(*HexLiteral[c4,Quotes->4]*)


(* ::Text:: *)
(*Error on the minimax approximation:*)


(* ::Input:: *)
(*End[]*)


(* ::Section:: *)
(*Error Analysis*)


(* ::Input:: *)
(*Options[mullerE]={UseFMA->False};*)
(*mullerE[\[Epsilon]1_Real,OptionsPattern[]]:=Module[{k=Floor[-Log2[\[Epsilon]1]-M],e},Assert[k>=3];e=If[OptionValue[UseFMA],1,(1-2^-M)^-1](1+2^(M+1) \[Epsilon]1/(1-\[Epsilon]1-2^(-k+1)));Assert[e<=2];e]*)


(* ::Subsection::Closed:: *)
(*Reduced Angle*)


(* ::Input:: *)
(*\[Zeta]0Interval=Interval[{-2^(-\[Kappa]3-M),2^(-\[Kappa]3-M)}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]0Interval]*)


(* ::Subsection::Closed:: *)
(*Sin Near Zero*)


(* ::Input:: *)
(*Begin["ErrorSinNearZero`"]*)


(* ::Input:: *)
(*x0Max=Max[x0Interval]*)


(* ::Text:: *)
(*Error on the minimax approximation:*)


(* ::Input:: *)
(*\[Eta]=Abs[sin0ApproximationResult[[2,2]]]*)


(* ::Text:: *)
(*Error due to the floating-point evaluation:*)


(* ::Input:: *)
(*\[Zeta]1=Interval[{-\[Eta],\[Eta]}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]1]*)


(* ::Input:: *)
(*\[Zeta]2=IEEEEvaluateWithRelativeError[sin0Polynomial[x0Interval^2]][[2]];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]2]*)


(* ::Input:: *)
(*\[Zeta]3=IEEEEvaluateWithRelativeError[x0Interval^3][[2]];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]3]*)


(* ::Input:: *)
(*\[Delta]1=uInterval;*)


(* ::Input:: *)
(*t0[x\:0303_]:=(1+\[Zeta]1)Sin[x\:0303]-x\:0303*)


(* ::Input:: *)
(*t1t2[x\:0303_]:=t0[x\:0303] (1+\[Zeta]2)(1+\[Zeta]3)*)


(* ::Input:: *)
(*t3[x\:0303_,\[Delta]x\:0303_]:=(t1t2[x\:0303]+\[Delta]x\:0303)(1+\[Delta]1)*)


(* ::Input:: *)
(*sin0ImplementationRelativeError[x\:0303_,\[Delta]x\:0303_]:=(x\:0303+t3[x\:0303,\[Delta]x\:0303])/Sin[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Text:: *)
(*The relative error is monotonic over the domain of interest:*)


(* ::Input:: *)
(*Plot3D[{Min[sin0ImplementationRelativeError[x,\[Delta]x]],Max[sin0ImplementationRelativeError[x,\[Delta]x]]},{x,0,x0Max},{\[Delta]x,-x0Max \[GothicU][1/2],x0Max \[GothicU][1/2]},RegionFunction->Function[{x,\[Delta]x},-x \[GothicU][1/2]<\[Delta]x<x \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*smol=1*^-100;*)


(* ::Input:: *)
(*corners={{smol,-smol \[GothicU][1/2]},{smol,smol \[GothicU][1/2]},{x0Max,-x0Max \[GothicU][1/2]},{x0Max,x0Max \[GothicU][1/2]}};*)


(* ::Input:: *)
(*relativeErrorsAtCorners=Block[{$MaxExtraPrecision=1000},Map[Apply[sin0ImplementationRelativeError,#]&,corners]];*)


(* ::Input:: *)
(*sin0ImplementationMaxRelativeError=Max[Abs/@relativeErrorsAtCorners]*)


(* ::Input:: *)
(*Log2[sin0ImplementationMaxRelativeError]*)


(* ::Text:: *)
(*This error preserves 17.6 bits after the mantissa of the result, which is sufficient for ensuring correct rounding with high probability.*)


(* ::Text:: *)
(*Rounding test:*)


(* ::Input:: *)
(*e=mullerE[sin0ImplementationMaxRelativeError,UseFMA->True]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Text:: *)
(*Proof that the term in \[Delta]x x^2 doesn't matter:*)


(* ::Input:: *)
(*Collect[Normal[Series[Sin[u],{u,0,5}]]/.u->(x+\[Delta]x),\[Delta]x]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*sin0NextTermRelativeSize[x\:0303_,\[Delta]x\:0303_]:=(-x\:0303^2\[Delta]x\:0303/2)/Sin[x\:0303+\[Delta]x\:0303]*)


(* ::Input:: *)
(*Plot3D[sin0NextTermRelativeSize[x,\[Delta]x],{x,0,x0Max},{\[Delta]x,-x0Max \[GothicU][1/2],x0Max \[GothicU][1/2]},RegionFunction->Function[{x,\[Delta]x},-x \[GothicU][1/2]<\[Delta]x<x \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}}]*)


(* ::Input:: *)
(*nextTermsAtCorners=Map[Apply[sin0NextTermRelativeSize,#]&,corners];*)


(* ::Input:: *)
(*sin0NextTermMaxRelativeSize=N[Max[Abs/@nextTermsAtCorners],20]*)


(* ::Input:: *)
(*Log2[sin0NextTermMaxRelativeSize]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsection::Closed:: *)
(*Around Table Entries*)


(* ::Text:: *)
(*These analyses  will need to be documented\:202f!*)


(* ::Text:: *)
(*\[Delta]x is less than half a ULP of \[Pi]/4:*)


(* ::Input:: *)
(*\[Delta]xInterval=Interval[{-\[GothicU][\[Pi]/4]/2,\[GothicU][\[Pi]/4]/2}];*)


(* ::Text:: *)
(*Error on the minimax approximations:*)


(* ::Input:: *)
(*\[Eta]s=Abs[sinApproximationResult[[2,2]]]*)


(* ::Input:: *)
(*\[Eta]c=Abs[cosApproximationResult[[2,2]]]*)


(* ::Input:: *)
(*\[Zeta]1=Interval[{-\[Eta]s,\[Eta]s}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]1]*)


(* ::Input:: *)
(*\[Zeta]2=Interval[{-\[Eta]c,\[Eta]c}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]2]*)


(* ::Text:: *)
(*Error due to the floating-point evaluation:*)


(* ::Input:: *)
(*sinApproximationIEEEIntervals=IEEEEvaluateWithRelativeError[sinPolynomial[hInterval^2]]*)


(* ::Input:: *)
(*cosApproximationIEEEIntervals=IEEEEvaluateWithRelativeError[ cosPolynomial[hInterval^2]]*)


(* ::Input:: *)
(*\[Zeta]3=sinApproximationIEEEIntervals[[2]];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]3]*)


(* ::Input:: *)
(*\[Zeta]4=cosApproximationIEEEIntervals[[2]];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]4]*)


(* ::Subsubsection::Closed:: *)
(*Sin*)


(* ::Input:: *)
(*Begin["ErrorSinAroundTableEntries`"]*)


(* ::Input:: *)
(*\[Delta]1=uInterval;*)


(* ::Input:: *)
(*\[Delta]2=uInterval;*)


(* ::Input:: *)
(*\[Delta]3=uInterval;*)


(* ::Input:: *)
(*\[Delta]4=uInterval;*)


(* ::Input:: *)
(*\[Delta]5=uInterval;*)


(* ::Input:: *)
(*\[Delta]6=uInterval;*)


(* ::Input:: *)
(*\[Delta]7=uInterval;*)


(* ::Input:: *)
(*\[Delta]8=uInterval;*)


(* ::Input:: *)
(*\[Delta]9=uInterval;*)


(* ::Input:: *)
(*\[Delta]10=uInterval;*)


(* ::Input:: *)
(*(*\[Delta]1=.;\[Delta]2=.;\[Delta]3=.;\[Delta]4=.;\[Delta]5=.;\[Delta]6=.;\[Delta]7=.;\[Delta]8=.;*)*)


(* ::Input:: *)
(*t0[h_,sk_,ck_]:=Hold[CorrectlyRound[ck]]h+Hold[CorrectlyRound[sk]]*)


(* ::Input:: *)
(*t1[h_]:=(1+\[Zeta]3)(1+\[Zeta]1)(Sin[h]-h)/h^3*)


(* ::Input:: *)
(*t2[h_]:=(1+\[Zeta]4)(1+\[Zeta]2)(Cos[h]-1)/h^2*)


(* ::Input:: *)
(*t3[h_,\[Delta]x\:0303_]:=h(2 \[Delta]x\:0303+h)(1+\[Delta]1)(1+\[Delta]2)*)


(* ::Input:: *)
(*t4[h_,\[Delta]x\:0303_]:=h^3(1+\[Delta]3)(1+\[Delta]4)*)


(* ::Input:: *)
(*t5[h_,\[Delta]x\:0303_,sk_]:=Hold[CorrectlyRound[sk]] t2[h]t3[h,\[Delta]x\:0303](1+\[Delta]5)(1+\[Delta]6)*)


(* ::Input:: *)
(*t6[h_,\[Delta]x\:0303_]:=t1[h]t4[h,\[Delta]x\:0303](1+\[Delta]7)*)


(* ::Input:: *)
(*t7[h_,\[Delta]x\:0303_,sk_,ck_]:=(Hold[CorrectlyRound[ck]]t6[h,\[Delta]x\:0303]+t5[h,\[Delta]x\:0303,sk])(1+\[Delta]8)*)


(* ::Input:: *)
(*t8[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=(Hold[CorrectlyRound[ck]] \[Delta]x\:0303+\[Epsilon] t0[h,sk,ck])(1+\[Delta]9)*)


(* ::Input:: *)
(*t9[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=(t7[h,\[Delta]x\:0303,sk,ck]+t8[h,\[Delta]x\:0303,sk,ck,\[Epsilon]])(1+\[Delta]10)*)


(* ::Input:: *)
(*t10[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=t0[h,sk,ck](1-\[Epsilon])+t9[h,\[Delta]x\:0303,sk,ck,\[Epsilon]]*)


(* ::Input:: *)
(*t11[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=ReleaseHold[CoefficientList[Collect[t10[h,\[Delta]x\:0303,sk,ck,\[Epsilon]],\[Epsilon]],\[Epsilon]]]*)


(* ::Input:: *)
(*t12[h_,\[Delta]x\:0303_,sk_,ck_]:=Module[{\[Epsilon],cl=t11[h,\[Delta]x\:0303,sk,ck,\[Epsilon]]},cl[[1]]+uInterval cl[[2]]]*)


(* ::Input:: *)
(*sinImplementationRelativeError[x\:0303_,\[Delta]x\:0303_,xk_,sk_,ck_]:=t12[x\:0303-xk,\[Delta]x\:0303,sk,ck]/Sin[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Input:: *)
(*m=Max[Abs[accurateTablesHIntervals[[65]]]]*)


(* ::Input:: *)
(*at=accurateTables[65]*)


(* ::Input:: *)
(*(*\[Zeta]0Interval=Interval[{0,0}]*)*)


(* ::Input:: *)
(*(*\[Delta]11=Interval[{0,0}];\[Delta]10=\[Delta]11;\[Delta]9=\[Delta]10;\[Delta]8=\[Delta]9;\[Delta]7=\[Delta]8;\[Delta]6=\[Delta]7;\[Delta]5=\[Delta]6;\[Delta]4=\[Delta]5;\[Delta]3=\[Delta]4;\[Delta]2=\[Delta]3;\[Delta]1=\[Delta]2;*)*)


(* ::Input:: *)
(*(*\[Zeta]4=Interval[{0,0}];\[Zeta]3=\[Zeta]4;\[Zeta]2=\[Zeta]3;\[Zeta]1=\[Zeta]2;*)*)


(* ::Input:: *)
(*sinImplementationRelativeError[at[[1]]+m,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]*)


(* ::Input:: *)
(*Plot[{Min[sinImplementationRelativeError[at[[1]]+h,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]],Max[sinImplementationRelativeError[at[[1]]+h,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]]},{h,0,m},WorkingPrecision->30]*)


(* ::Input:: *)
(*Plot3D[{Min[sinImplementationRelativeError[at[[1]]+h,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]],Max[sinImplementationRelativeError[at[[1]]+h,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]]},{h,0,m},{\[Delta]x\:0303,-m \[GothicU][1/2],m \[GothicU][1/2]},WorkingPrecision->30,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*smol=1*^-100;*)


(* ::Input:: *)
(*sinImplementationMaxRelativeErrorPerInterval[i_]:=*)
(*Module[*)
(*{at=accurateTables[i],m=Max[Abs[accurateTablesHIntervals[[i]]]],corners,r},*)
(*corners={{smol,-m \[GothicU][1/2]},{smol,m \[GothicU][1/2]},{m,-m \[GothicU][1/2]},{m,m \[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[sinImplementationRelativeError[at[[1]]+#[[1]],#[[2]],at[[1]],at[[2]],at[[3]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{i,sinImplementationMaxRelativeErrorPerInterval[i]},{i,1,accurateTablesMaxIndex}]]*)


(* ::Text:: *)
(*The dispersion of the error is largely due to the bit pattern of sk after the 18th accurate bit:*)


(* ::Input:: *)
(*Bits[accurateTables[65][[2]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[89][[2]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[283][[2]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[399][[2]],30]*)


(* ::Input:: *)
(*sinImplementationMaxRelativeError=Max[Table[sinImplementationMaxRelativeErrorPerInterval[i],{i,1,accurateTablesMaxIndex}]];*)


(* ::Input:: *)
(*Log2[sinImplementationMaxRelativeError]*)


(* ::Text:: *)
(*Rounding test:*)


(* ::Input:: *)
(*e=mullerE[sinImplementationMaxRelativeError,UseFMA->True];*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Text:: *)
(*Proof that the terms in h^2\[Delta]h and above can be ignored:*)


(* ::Input:: *)
(*su=Collect[Normal[Series[Sin[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*cu=Collect[Normal[Series[Cos[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*Collect[TrigExpand[Sin[xk+u]]/.{Sin[u]->su,Cos[u]->cu},{\[Delta]x,h}]*)


(* ::Input:: *)
(*sinNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=(\[Delta]x\:0303 (x\:0303-xk)^2Cos[xk]/2)/Sin[x\:0303+\[Delta]x\:0303]*)


(* ::Input:: *)
(*Plot3D[sinNextTermRelativeSize[at[[1]]+h,\[Delta]x\:0303,at[[1]]],{h,0,m},{\[Delta]x\:0303,-m \[GothicU][1/2],m \[GothicU][1/2]},WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotRange->Full]*)


(* ::Input:: *)
(*nextTermsAtCorners[i_]:=Module[*)
(*{at=accurateTables[i],m=Max[Abs[accurateTablesHIntervals[[i]]]],corners,r},*)
(*corners={{smol,-m \[GothicU][1/2]},{smol,m \[GothicU][1/2]},{m,-m \[GothicU][1/2]},{m,m \[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[sinNextTermRelativeSize[at[[1]]+#[[1]],#[[2]],at[[1]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{i,nextTermsAtCorners[i]},{i,1,accurateTablesMaxIndex}]]*)


(* ::Input:: *)
(*sinNextTermMaxRelativeSize=Max[Table[nextTermsAtCorners[i],{i,1,accurateTablesMaxIndex}]];*)


(* ::Input:: *)
(*N[Log2[sinNextTermMaxRelativeSize]]*)


(* ::Input:: *)
(*End[]*)


(* ::Subsubsection::Closed:: *)
(*Cos*)


(* ::Input:: *)
(*Begin["ErrorsCosAroundTableEntries`"]*)


(* ::Input:: *)
(*\[Delta]1=uInterval;*)


(* ::Input:: *)
(*\[Delta]2=uInterval;*)


(* ::Input:: *)
(*\[Delta]3=uInterval;*)


(* ::Input:: *)
(*\[Delta]4=uInterval;*)


(* ::Input:: *)
(*\[Delta]5=uInterval;*)


(* ::Input:: *)
(*\[Delta]6=uInterval;*)


(* ::Input:: *)
(*\[Delta]7=uInterval;*)


(* ::Input:: *)
(*\[Delta]8=uInterval;*)


(* ::Input:: *)
(*\[Delta]9=uInterval;*)


(* ::Input:: *)
(*\[Delta]10=uInterval;*)


(* ::Input:: *)
(*(*\[Delta]1=.;\[Delta]2=.;\[Delta]3=.;\[Delta]4=.;\[Delta]5=.;\[Delta]6=.;\[Delta]7=.;\[Delta]8=.;*)*)


(* ::Input:: *)
(*t0[h_,sk_,ck_]:=-Hold[CorrectlyRound[sk]]h+Hold[CorrectlyRound[ck]]*)


(* ::Input:: *)
(*t1[h_]:=(1+\[Zeta]3)(1+\[Zeta]1)(Sin[h]-h)/h^3*)


(* ::Input:: *)
(*t2[h_]:=(1+\[Zeta]4)(1+\[Zeta]2)(Cos[h]-1)/h^2*)


(* ::Input:: *)
(*t3[h_,\[Delta]x\:0303_]:=h(2 \[Delta]x\:0303+h)(1+\[Delta]1)(1+\[Delta]2)*)


(* ::Input:: *)
(*t4[h_,\[Delta]x\:0303_]:=h^3(1+\[Delta]3)(1+\[Delta]4)*)


(* ::Input:: *)
(*t5[h_,\[Delta]x\:0303_,ck_]:=Hold[CorrectlyRound[ck]] t3[h,\[Delta]x\:0303]t2[h](1+\[Delta]5)(1+\[Delta]6)*)


(* ::Input:: *)
(*t6[h_,\[Delta]x\:0303_]:=t4[h,\[Delta]x\:0303]t1[h](1+\[Delta]7)*)


(* ::Input:: *)
(*t7[h_,\[Delta]x\:0303_,sk_,ck_]:=(-Hold[CorrectlyRound[sk]]t6[h,\[Delta]x\:0303]+t5[h,\[Delta]x\:0303,ck])(1+\[Delta]8)*)


(* ::Input:: *)
(*t8[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=(-Hold[CorrectlyRound[sk]] \[Delta]x\:0303+\[Epsilon] t0[h,sk,ck])(1+\[Delta]9)*)


(* ::Input:: *)
(*t9[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=(t7[h,\[Delta]x\:0303,sk,ck]+t8[h,\[Delta]x\:0303,sk,ck,\[Epsilon]])(1+\[Delta]10)*)


(* ::Input:: *)
(*t10[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=t0[h,sk,ck](1-\[Epsilon])+t9[h,\[Delta]x\:0303,sk,ck,\[Epsilon]]*)


(* ::Input:: *)
(*t11[h_,\[Delta]x\:0303_,sk_,ck_,\[Epsilon]_]:=ReleaseHold[CoefficientList[Collect[t10[h,\[Delta]x\:0303,sk,ck,\[Epsilon]],\[Epsilon]],\[Epsilon]]]*)


(* ::Input:: *)
(*t12[h_,\[Delta]x\:0303_,sk_,ck_]:=Module[{\[Epsilon],cl=t11[h,\[Delta]x\:0303,sk,ck,\[Epsilon]]},cl[[1]]+uInterval cl[[2]]]*)


(* ::Input:: *)
(*cosImplementationRelativeError[x\:0303_,\[Delta]x\:0303_,xk_,sk_,ck_]:=t12[x\:0303-xk,\[Delta]x\:0303,sk,ck]/Cos[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Input:: *)
(*m=Max[Abs[accurateTablesHIntervals[[65]]]]*)


(* ::Input:: *)
(*at=accurateTables[396]*)


(* ::Input:: *)
(*(*\[Zeta]0Interval=Interval[{0,0}]*)*)


(* ::Input:: *)
(*(*\[Delta]11=Interval[{0,0}];\[Delta]10=\[Delta]11;\[Delta]9=\[Delta]10;\[Delta]8=\[Delta]9;\[Delta]7=\[Delta]8;\[Delta]6=\[Delta]7;\[Delta]5=\[Delta]6;\[Delta]4=\[Delta]5;\[Delta]3=\[Delta]4;\[Delta]2=\[Delta]3;\[Delta]1=\[Delta]2;*)*)


(* ::Input:: *)
(*(*\[Zeta]4=Interval[{0,0}];\[Zeta]3=\[Zeta]4;\[Zeta]2=\[Zeta]3;\[Zeta]1=\[Zeta]2;*)*)


(* ::Input:: *)
(*cosImplementationRelativeError[at[[1]]+m,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]*)


(* ::Input:: *)
(*Plot[{Min[cosImplementationRelativeError[at[[1]]+h,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]],Max[cosImplementationRelativeError[at[[1]]+h,m \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]]},{h,0,m},WorkingPrecision->30]*)


(* ::Input:: *)
(*Plot3D[{Min[cosImplementationRelativeError[at[[1]]+h,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]],Max[cosImplementationRelativeError[at[[1]]+h,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]]},{h,0,m},{\[Delta]x\:0303,-m \[GothicU][1/2],m \[GothicU][1/2]},WorkingPrecision->30,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*smol=1*^-100;*)


(* ::Input:: *)
(*cosImplementationMaxRelativeErrorPerInterval[i_]:=*)
(*Module[*)
(*{at=accurateTables[i],m=Max[Abs[accurateTablesHIntervals[[i]]]],corners,r},*)
(*corners={{smol,-m \[GothicU][1/2]},{smol,m \[GothicU][1/2]},{m,-m \[GothicU][1/2]},{m,m \[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[cosImplementationRelativeError[at[[1]]+#[[1]],#[[2]],at[[1]],at[[2]],at[[3]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{i,cosImplementationMaxRelativeErrorPerInterval[i]},{i,1,accurateTablesMaxIndex}]]*)


(* ::Text:: *)
(*The dispersion of the error is largely due to the bit pattern of ck after the 18th accurate bit:*)


(* ::Input:: *)
(*Bits[accurateTables[387][[3]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[13][[3]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[396][[3]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[83][[3]],30]*)


(* ::Input:: *)
(*cosImplementationMaxRelativeError=Max[Table[cosImplementationMaxRelativeErrorPerInterval[i],{i,1,accurateTablesMaxIndex}]];*)


(* ::Input:: *)
(*Log2[cosImplementationMaxRelativeError]*)


(* ::Text:: *)
(*Rounding test:*)


(* ::Input:: *)
(*e=mullerE[cosImplementationMaxRelativeError,UseFMA->True];*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Text:: *)
(*Proof that the terms in h^2 \[Delta]h and above can be ignored:*)


(* ::Input:: *)
(*su=Collect[Normal[Series[Sin[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*cu=Collect[Normal[Series[Cos[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*Collect[TrigExpand[Cos[xk+u]]/.{Sin[u]->su,Cos[u]->cu},{\[Delta]x,h}]*)


(* ::Input:: *)
(*cosNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=(\[Delta]x\:0303 (x\:0303-xk)^2Sin[xk]/2)/Cos[x\:0303+\[Delta]x\:0303]*)


(* ::Input:: *)
(*Plot3D[cosNextTermRelativeSize[at[[1]]+h,\[Delta]x\:0303,at[[1]]],{h,0,m},{\[Delta]x\:0303,-m \[GothicU][1/2],m \[GothicU][1/2]},WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotRange->Full]*)


(* ::Input:: *)
(*nextTermsAtCorners[i_]:=Module[*)
(*{at=accurateTables[i],m=Max[Abs[accurateTablesHIntervals[[i]]]],corners,r},*)
(*corners={{smol,-m \[GothicU][1/2]},{smol,m \[GothicU][1/2]},{m,-m \[GothicU][1/2]},{m,m \[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[cosNextTermRelativeSize[at[[1]]+#[[1]],#[[2]],at[[1]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{i,nextTermsAtCorners[i]},{i,1,accurateTablesMaxIndex}]]*)


(* ::Input:: *)
(*cosNextTermMaxRelativeSize=Max[Table[nextTermsAtCorners[i],{i,1,accurateTablesMaxIndex}]];*)


(* ::Input:: *)
(*N[Log2[cosNextTermMaxRelativeSize]]*)


(* ::Input:: *)
(*End[]*)
