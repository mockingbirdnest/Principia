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
(*Find the bounds over which we'll have to evaluate the polynomial, based on the excursion away from multiples of 1/512, except for k=1:*)


(* ::Input:: *)
(*accurateTablesXIntervals=Table[Interval[{If[k==1,(1+1/8)accurateTablesStep/2,(2k-1)accurateTablesStep/2],(2k+1)accurateTablesStep/2}],{k,1,accurateTablesMaxIndex}];*)


(* ::Input:: *)
(*accurateTablesHIntervals=Table[accurateTablesXIntervals[[k]]-accurateTables[k][[1]],{k,1,accurateTablesMaxIndex}];*)


(* ::Input:: *)
(*{hLB,hUB}={Min[Min/@accurateTablesHIntervals],Max[Max/@accurateTablesHIntervals]};*)


(* ::Input:: *)
(*hInterval=Interval[{hLB,hUB}];*)


(* ::Input:: *)
(*N[Log2[Abs[{hLB,hUB}]]]*)


(* ::Input:: *)
(* N[Log2[Abs[{hLB,hUB}]-accurateTablesStep/2]]*)


(* ::Text:: *)
(*Perturbations:*)


(* ::Input:: *)
(*accurateTablesPerturbations=Table[Representation[accurateTables[i][[1]]]-Representation[i/512],{i,1,accurateTablesMaxIndex}];*)


(* ::Input:: *)
(*accurateTablesPerturbedBits=Map[N[Log2[Abs[#]]]&,accurateTablesPerturbations];*)


(* ::Input:: *)
(*PositionLargest[accurateTablesPerturbedBits]*)


(* ::Input:: *)
(*PositionLargest[accurateTablesPerturbedBits[[2;;]]]*)


(* ::Input:: *)
(*PositionSmallest[accurateTablesPerturbedBits[[2;;]]]*)


(* ::Input:: *)
(*accurateTablesPerturbedBits[[1]]*)


(* ::Input:: *)
(*accurateTablesPerturbedBits[[8]]*)


(* ::Input:: *)
(*accurateTablesPerturbedBits[[57]]*)


(* ::Text:: *)
(*Check the Sterbenz condition for computing sk+ck h exactly (the subtraction h' - sk is exact):*)


(* ::Input:: *)
(*AllTrue[Table[Module[{t=accurateTables[k],sk,ck,h,h\[Prime]},*)
(*sk=t[[2]];ck=t[[3]];h=accurateTablesHIntervals[[k]];h\[Prime]=(sk+ck h)(1+uInterval);sk/2<=h\[Prime]<=2sk],{k,1,accurateTablesMaxIndex}],TrueQ]*)


(* ::Text:: *)
(*Similarly the Sterbenz condition for computing ck - h sk exactly (the subtraction h'- ck is exact):*)


(* ::Input:: *)
(*AllTrue[Table[Module[{t=accurateTables[k],sk,ck,h,h\[Prime]},*)
(*sk=t[[2]];ck=t[[3]];h=accurateTablesHIntervals[[k]];h\[Prime]=(ck-h sk)(1+uInterval);ck/2<=h\[Prime]<=2 ck],{k,1,accurateTablesMaxIndex}],TrueQ]*)


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


(* ::Input:: *)
(*N[Log2[{Abs[\[Pi]/2-C1-\[Delta]C1],2^(\[Kappa]\[Prime]1-M-1)\[GothicU][\[Pi]/2]}],20]*)


(* ::Text:: *)
(*That's probably because there are three zeroes after \[Delta]C1:*)


(* ::Input:: *)
(*Bits[\[Pi]/2-C1]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C1]<2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C1],2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*N[Log2[{Abs[\[Delta]C1],2^\[Kappa]\[Prime]1(1+2^(-M-1))\[GothicU][\[Pi]/2]}],20]*)


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


(* ::Input:: *)
(*N[Log2[{Abs[\[Pi]/2-C2-C\[Prime]2-\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-2 M-1)\[GothicU][\[Pi]/2]}],20]*)


(* ::Text:: *)
(*That's because the first bit after \[Delta]C2 is one:*)


(* ::Input:: *)
(*Bits[\[Pi]/2-C2-C\[Prime]2]*)


(* ::Input:: *)
(*Assert[Abs[\[Delta]C2]<2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M)(1+2^(-M-1))\[GothicU][\[Pi]/2]]*)


(* ::Input:: *)
(*N[{Abs[\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M)(1+2^(-M-1))\[GothicU][\[Pi]/2]},20]*)


(* ::Input:: *)
(*N[Log2[{Abs[\[Delta]C2],2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M)(1+2^(-M-1))\[GothicU][\[Pi]/2]}],20]*)


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


(* ::Input:: *)
(*N[Log2[{Max[reductionInterval],2^(\[Kappa]1+\[Kappa]\[Prime]1-M+1)\[GothicU][\[Pi]/2]}],20]*)


(* ::Text:: *)
(*Threshold for dangerous rounding:*)


(* ::Input:: *)
(*angleReducedThreshold = 2^(\[Kappa]1+\[Kappa]\[Prime]1+\[Kappa]3-M+2)*)


(* ::Input:: *)
(*Denominator[angleReducedThreshold]//Log2*)


(* ::Input:: *)
(*N[angleReducedThreshold]*)


(* ::Text:: *)
(*Tighter threshold possible:*)


(* ::Input:: *)
(*N[2^(\[Kappa]3+M)Max[reductionInterval]]*)


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
(*\[Zeta]2=2^(1-2M);*)


(* ::Text:: *)
(*Error on the overall reduction:*)


(* ::Input:: *)
(*reductionInterval=(\[Zeta]1 nInterval + \[Delta]yInterval[[2]])(1+\[Zeta]2)+Interval[{-\[Pi]/4-misroundingBound,\[Pi]/4+misroundingBound}]\[Zeta]2;*)


(* ::Input:: *)
(*Assert[Max[reductionInterval]<2^(\[Kappa]2-2M)(2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2+1)+5)\[GothicU][\[Pi]/2]+2^(-2M)(\[Pi]/2)]*)


(* ::Input:: *)
(*N[{Max[reductionInterval],2^(\[Kappa]2-2M)(2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2+1)+5)\[GothicU][\[Pi]/2]+2^(-2M)(\[Pi]/2)},20]*)


(* ::Input:: *)
(*N[Log2[{Max[reductionInterval],2^(\[Kappa]2-2M)(2^(\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2+1)+5)\[GothicU][\[Pi]/2]+2^(-2M)(\[Pi]/2)}],20]*)


(* ::Text:: *)
(*Threshold for dangerous rounding:*)


(* ::Input:: *)
(*angleReducedThreshold = 2^(\[Kappa]3-M)(2^(\[Kappa]2+\[Kappa]\[Prime]2+\[Kappa]\[DoublePrime]2-M+2)+2)*)


(* ::Input:: *)
(*Denominator[angleReducedThreshold]//Log2*)


(* ::Input:: *)
(*N[angleReducedThreshold]*)


(* ::Text:: *)
(*Tighter threshold possible:*)


(* ::Input:: *)
(*N[2^(\[Kappa]3+M)Max[reductionInterval]]*)


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
(*x0Max=Min[accurateTablesXIntervals[[1]]];*)


(* ::Input:: *)
(*x0Min=-x0Max;*)


(* ::Input:: *)
(*x0Interval=Interval[{x0Min,x0Max}]*)


(* ::Text:: *)
(*The value of the error function at 0 must be chosen so that we have proper convergence (which we can check by looking for equioscillation in the graph).*)


(* ::Input:: *)
(*sin0ApproximationResult=GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1*^6,Sin[t]/t^3]},*)
(*{t,{0,x0Max},1,0},*)
(*x,WorkingPrecision->30]*)


(* ::Input:: *)
(*sin0Polynomial=Function[u, Evaluate[HornerForm[sin0ApproximationResult[[2,1]]/.x->u]]]*)


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


(* ::Text:: *)
(*The value of the error function at 0 must be chosen so that we have proper convergence (which we can check by looking for equioscillation in the graph).*)


(* ::Input:: *)
(*sinApproximationResult=GeneralMiniMaxApproximation[*)
(*{t^2,If[t==0,-1/6,sinFn[t]],If[t==0,1*^6,Sin[t]/t^3]},*)
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
(*cosPolynomial=Function[u, Evaluate[HornerForm[ cosApproximationResult[[2,1]]/.x->u]]]*)


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


(* ::Input:: *)
(*End[]*)


(* ::Section:: *)
(*Error Analysis*)


(* ::Text:: *)
(*This section analyses the relative error of the various algorithms and computes Muller's e factor.*)


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
(*\[Xi]=Abs[sin0ApproximationResult[[2,2]]]*)


(* ::Text:: *)
(*Error due to the floating-point evaluation:*)


(* ::Input:: *)
(*\[Zeta]1=Interval[{-\[Xi],\[Xi]}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]1]*)


(* ::Input:: *)
(*\[Zeta]2=IEEEEvaluateWithRelativeError[sin0Polynomial[x0Interval^2]][[2]];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]2]*)


(* ::Input:: *)
(*\[Delta]1=uInterval;*)


(* ::Input:: *)
(*\[Delta]2=uInterval;*)


(* ::Input:: *)
(*\[Delta]3=uInterval;*)


(* ::Input:: *)
(*\[Delta]4=uInterval;*)


(* ::Input:: *)
(*t1[x\:0303_]:=(1+\[Zeta]2)((1+\[Zeta]1)Sin[x\:0303]-x\:0303)/x\:0303^3*)


(* ::Input:: *)
(*t2[x\:0303_]:=x\:0303^3(1+\[Delta]1)(1+\[Delta]2)*)


(* ::Input:: *)
(*t3[x\:0303_,\[Delta]x\:0303_]:=(t1[x\:0303]t2[x\:0303](1+\[Delta]3)+\[Delta]x\:0303)(1+\[Delta]4)*)


(* ::Input:: *)
(*sin0ImplementationRelativeError[x\:0303_,\[Delta]x\:0303_]:=(x\:0303+t3[x\:0303,\[Delta]x\:0303])/Sin[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Input:: *)
(*Plot[{Min[sin0ImplementationRelativeError[x\:0303,0]],Max[sin0ImplementationRelativeError[x\:0303,0]]},{x\:0303,0,x0Max},WorkingPrecision->30]*)


(* ::Text:: *)
(*The relative error is monotonic over the domain of interest:*)


(* ::Input:: *)
(*smol=1*^-100;*)


(* ::Input:: *)
(*Plot3D[{Min[sin0ImplementationRelativeError[x\:0303,\[Delta]x\:0303]],Max[sin0ImplementationRelativeError[x\:0303,\[Delta]x\:0303]]},{x\:0303,smol,x0Max},{\[Delta]x\:0303,-x0Max \[GothicU][1/2],x0Max \[GothicU][1/2]},RegionFunction->Function[{x\:0303,\[Delta]x\:0303},-x\:0303 \[GothicU][1/2]<\[Delta]x\:0303<x\:0303 \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*corners={{smol,-smol \[GothicU][1/2]},{smol,smol \[GothicU][1/2]},{x0Max,-x0Max \[GothicU][1/2]},{x0Max,x0Max \[GothicU][1/2]}};*)


(* ::Input:: *)
(*relativeErrorsAtCorners=Block[{$MaxExtraPrecision=1000},Map[Apply[sin0ImplementationRelativeError,#]&,corners]];*)


(* ::Input:: *)
(*sin0ImplementationMaxRelativeError=Max[Abs/@relativeErrorsAtCorners]*)


(* ::Input:: *)
(*Log2[sin0ImplementationMaxRelativeError]*)


(* ::Text:: *)
(*Rounding test:*)


(* ::Input:: *)
(*e=mullerE[sin0ImplementationMaxRelativeError,UseFMA->False]*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Subsubsection::Closed:: *)
(*Dominant Term of the Error*)


(* ::Input:: *)
(*errorExpression=Block[{\[Delta]1=d1,\[Delta]2=d2,\[Delta]3=d3,\[Delta]4=d4,\[Zeta]1=z1,\[Zeta]2=z2},t3[x\:0303,\[Delta]x\:0303]]*)


(* ::Input:: *)
(*vars={d1,d2,d3,d4,z1,z2,\[Delta]x\:0303};*)


(* ::Input:: *)
(*errorExpression1stOrder=Module[{alt=Alternatives@@vars},Collect[Expand[errorExpression]/.Times->times/.times[___,alt,___,alt,___]->0/.times->Times,vars]]*)


(* ::Subsubsection::Closed:: *)
(*Proof That the Term in \[Delta]x x^2 Does Not Matter*)


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
(*Error on the minimax approximations:*)


(* ::Input:: *)
(*\[Xi]s=Abs[sinApproximationResult[[2,2]]]*)


(* ::Input:: *)
(*Log2[\[Xi]s]*)


(* ::Input:: *)
(*\[Xi]c=Abs[cosApproximationResult[[2,2]]]*)


(* ::Input:: *)
(*Log2[\[Xi]c]*)


(* ::Input:: *)
(*\[Zeta]1=Interval[{-\[Xi]s,\[Xi]s}];*)


(* ::Input:: *)
(*binaryBounds[\[Zeta]1]*)


(* ::Input:: *)
(*\[Zeta]2=Interval[{-\[Xi]c,\[Xi]c}];*)


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


(* ::Input:: *)
(*\[Eta]=Interval[{-2^(-2M+1),2^(-2M+1)}]*)


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
(*\[Delta]11=uInterval;*)


(* ::Input:: *)
(*\[Delta]12=uInterval;*)


(* ::Input:: *)
(*t0[h_,sk_,ck_]:=(Hold[CorrectlyRound[ck]]h+Hold[CorrectlyRound[sk]])(1+\[Eta])*)


(* ::Input:: *)
(*t1[h_]:=((Sin[h](1+\[Zeta]1)-h)/h^3)(1+\[Zeta]3)*)


(* ::Input:: *)
(*t2[h_]:=((Cos[h]-1)/h^2)(1+\[Zeta]2)(1+\[Zeta]4)*)


(* ::Input:: *)
(*t3[h_,\[Delta]x\:0303_]:=h(2 \[Delta]x\:0303+h)(1+\[Delta]1)(1+\[Delta]2)*)


(* ::Input:: *)
(*t4[h_,\[Delta]x\:0303_]:=h^3(1+\[Delta]3)(1+\[Delta]4)*)


(* ::Input:: *)
(*t5[h_,\[Delta]x\:0303_,sk_]:=Hold[CorrectlyRound[sk]]t3[h,\[Delta]x\:0303] t2[h](1+\[Delta]5)(1+\[Delta]6)*)


(* ::Input:: *)
(*t6[h_,\[Delta]x\:0303_]:=(t4[h,\[Delta]x\:0303]t1[h](1+\[Delta]7)+\[Delta]x\:0303)(1+\[Delta]8)*)


(* ::Input:: *)
(*t7[h_,\[Delta]x\:0303_,sk_,ck_]:=(Hold[CorrectlyRound[ck]]t6[h,\[Delta]x\:0303](1+\[Delta]9)+t5[h,\[Delta]x\:0303,sk])(1+\[Delta]10)*)


(* ::Input:: *)
(*t8[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=(\[Delta]0 t0[h,sk,ck](1+\[Delta]11)+t7[h,\[Delta]x\:0303,sk,ck])(1+\[Delta]12)*)


(* ::Input:: *)
(*t9[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=t0[h,sk,ck](1-\[Delta]0)+t8[h,\[Delta]x\:0303,sk,ck,\[Delta]0]*)


(* ::Input:: *)
(*t10[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=ReleaseHold[CoefficientList[Collect[t9[h,\[Delta]x\:0303,sk,ck,\[Delta]0],\[Delta]0],\[Delta]0]]*)


(* ::Input:: *)
(*t11[h_,\[Delta]x\:0303_,sk_,ck_]:=Module[{\[Delta]0,cl=t10[h,\[Delta]x\:0303,sk,ck,\[Delta]0]},cl[[1]]+uInterval cl[[2]]]*)


(* ::Input:: *)
(*sinImplementationRelativeError[x\:0303_,\[Delta]x\:0303_,xk_,sk_,ck_]:=t11[x\:0303-xk,\[Delta]x\:0303,sk,ck]/Sin[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Input:: *)
(*Block[{\[Eta]=eta,\[Delta]1=d1,\[Delta]2=d2,\[Delta]3=d3,\[Delta]4=d4,\[Delta]5=d5,\[Delta]6=d6,\[Delta]7=d7,\[Delta]8=d8,\[Delta]9=d9,\[Delta]10=d10,\[Delta]11=d11,\[Delta]12=d12},CoefficientList[Collect[t9[h,\[Delta]x\:0303,sk,ck,\[Delta]0],\[Delta]0],\[Delta]0][[2]]]*)


(* ::Input:: *)
(*xi=accurateTablesXIntervals[[1]];*)


(* ::Input:: *)
(*at=accurateTables[1];*)


(* ::Input:: *)
(*Plot[{Min[sinImplementationRelativeError[x\:0303,Max[xi] \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]],Max[sinImplementationRelativeError[x\:0303,Max[xi] \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]]},{x\:0303,Min[xi],Max[xi]},WorkingPrecision->30,PlotRange->Full]*)


(* ::Input:: *)
(*Plot3D[{Min[sinImplementationRelativeError[x\:0303,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]],Max[sinImplementationRelativeError[x\:0303,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]]},{x\:0303,Min[xi],Max[xi]},{\[Delta]x\:0303,-Min[xi]\[GothicU][1/2],Max[xi] \[GothicU][1/2]},RegionFunction->Function[{x\:0303,\[Delta]x\:0303},-x\:0303 \[GothicU][1/2]<\[Delta]x\:0303<x\:0303 \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*sinImplementationMaxRelativeErrorPerInterval[k_]:=*)
(*Module[*)
(*{at=accurateTables[k],xi=accurateTablesXIntervals[[k]],corners,r},*)
(*corners={{Min[xi],-Min[xi]\[GothicU][1/2]},{Min[xi],Min[xi] \[GothicU][1/2]},{Max[xi],-Max[xi]\[GothicU][1/2]},{Max[xi],Max[xi]\[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[sinImplementationRelativeError[#[[1]],#[[2]],at[[1]],at[[2]],at[[3]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{k,sinImplementationMaxRelativeErrorPerInterval[k]},{k,1,accurateTablesMaxIndex}]]*)


(* ::Text:: *)
(*The dispersion of the error is largely due to the bit pattern of sk after the 18th accurate bit:*)


(* ::Input:: *)
(*Bits[accurateTables[17][[2]],30]*)


(* ::Input:: *)
(*Bits[accurateTables[65][[2]],30]*)


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
(*e=mullerE[sinImplementationMaxRelativeError,UseFMA->False];*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Subsubsubsection::Closed:: *)
(*Dominant Term of the Error*)


(* ::Input:: *)
(*errorExpression=Block[{\[Eta]=eta,\[Delta]1=d1,\[Delta]2=d2,\[Delta]3=d3,\[Delta]4=d4,\[Delta]5=d5,\[Delta]6=d6,\[Delta]7=d7,\[Delta]8=d8,\[Delta]9=d9,\[Delta]10=d10,\[Delta]11=d11,\[Delta]12=d12,\[Zeta]1=z1,\[Zeta]2=z2,\[Zeta]3=z3,\[Zeta]4=z4},CoefficientList[Collect[t9[h,\[Delta]x\:0303,sk,ck,\[Delta]0],\[Delta]0],\[Delta]0][[1]]/.Hold[CorrectlyRound[x_]]->x]*)


(* ::Input:: *)
(*vars={eta,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,z1,z2,z3,z4,\[Delta]0,\[Delta]x\:0303}*)


(* ::Input:: *)
(*errorExpression1stOrder=Module[{alt=Alternatives@@vars},Collect[Expand[errorExpression]/.Times->times/.times[___,alt,___,alt,___]->0/.times->Times,vars]]*)


(* ::Text:: *)
(*This is the main term:*)


(* ::Input:: *)
(*errorExpression1stOrder[[{2,9,10}]]*)


(* ::Input:: *)
(*xie=accurateTablesXIntervals[[65]];*)


(* ::Input:: *)
(*ate=accurateTables[65];*)


(* ::Input:: *)
(*errorExpression1stOrderInterval=(List@@errorExpression1stOrder)/.h->x\:0303-xk/.{d1->uInterval,d2->uInterval,d3->uInterval,d4->uInterval,d5->uInterval,d6->uInterval,d7->uInterval,d8->uInterval,d9->uInterval,d10->uInterval,d11->uInterval,d12->uInterval,eta->2^(-2M),z1->\[Zeta]1,z2->\[Zeta]2,z3->\[Zeta]3,z4->\[Zeta]4,\[Delta]x\:0303-> Max[xie]\[GothicU][1/2],xk->ate[[1]],sk->ate[[2]],ck->ate[[3]]};*)


(* ::Input:: *)
(*Plot[MinMax[(Plus@@errorExpression1stOrderInterval)/Sin[x\:0303+Max[xie]\[GothicU][1/2]+\[Zeta]0Interval x\:0303]-1],{x\:0303,Min[xie],Max[xie]},PlotRange->Full,WorkingPrecision->30,PlotLegends->Automatic]*)


(* ::Input:: *)
(*Table[Plot[MinMax[errorExpression1stOrderInterval[[i]]/Sin[x\:0303+Max[xie]\[GothicU][1/2]+\[Zeta]0Interval x\:0303]],{x\:0303,Min[xie],Max[xie]},PlotRange->Full,WorkingPrecision->30,PlotLegends->Automatic],{i,1,Length[errorExpression1stOrderInterval]}]*)


(* ::Text:: *)
(*The dominant term of the error is:*)


(* ::Input:: *)
(*errorExpression1stOrder[[7]]*)


(* ::Subsubsubsection::Closed:: *)
(*Proof That the Terms in h^2 \[Delta]x\:0303 and above Can Be Ignored*)


(* ::Input:: *)
(*su=Collect[Normal[Series[Sin[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*cu=Collect[Normal[Series[Cos[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*Collect[TrigExpand[Sin[xk+u]]/.{Sin[u]->su,Cos[u]->cu},{\[Delta]x,h}]*)


(* ::Input:: *)
(*sinNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=(\[Delta]x\:0303 (x\:0303-xk)^2Cos[xk]/2)/Sin[x\:0303+\[Delta]x\:0303]*)


(* ::Input:: *)
(*(*sinNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=-(\[Delta]x\:0303 (x\:0303-xk)Sin[xk])/Sin[x\:0303+\[Delta]x\:0303]*)*)


(* ::Input:: *)
(*Plot3D[sinNextTermRelativeSize[x\:0303,\[Delta]x\:0303,at[[1]]],{x\:0303,Min[xi],Max[xi]},{\[Delta]x\:0303,-Max[xi]\[GothicU][1/2],Max[xi] \[GothicU][1/2]},RegionFunction->Function[{x\:0303,\[Delta]x\:0303},-x\:0303 \[GothicU][1/2]<\[Delta]x\:0303<x\:0303 \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotRange->Full]*)


(* ::Input:: *)
(*nextTermsAtCorners[k_]:=Module[*)
(*{at=accurateTables[k],xi=accurateTablesXIntervals[[k]],corners,r},*)
(*corners={{Min[xi],-Min[xi]\[GothicU][1/2]},{Min[xi],Min[xi] \[GothicU][1/2]},{Max[xi],-Max[xi]\[GothicU][1/2]},{Max[xi],Max[xi]\[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[sinNextTermRelativeSize[#[[1]],#[[2]],at[[1]]]&,corners]];*)
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
(*\[Delta]11=uInterval;*)


(* ::Input:: *)
(*\[Delta]12=uInterval;*)


(* ::Input:: *)
(*t0[h_,sk_,ck_]:=(-Hold[CorrectlyRound[sk]]h+Hold[CorrectlyRound[ck]])(1+\[Eta])*)


(* ::Input:: *)
(*t1[h_]:=((Sin[h](1+\[Zeta]1)-h)/h^3)(1+\[Zeta]3)*)


(* ::Input:: *)
(*t2[h_]:=((Cos[h]-1)/h^2)(1+\[Zeta]2)(1+\[Zeta]4)*)


(* ::Input:: *)
(*t3[h_,\[Delta]x\:0303_]:=h(2 \[Delta]x\:0303+h)(1+\[Delta]1)(1+\[Delta]2)*)


(* ::Input:: *)
(*t4[h_,\[Delta]x\:0303_]:=h^3(1+\[Delta]3)(1+\[Delta]4)*)


(* ::Input:: *)
(*t5[h_,\[Delta]x\:0303_,ck_]:=Hold[CorrectlyRound[ck]] t3[h,\[Delta]x\:0303]t2[h](1+\[Delta]5)(1+\[Delta]6)*)


(* ::Input:: *)
(*t6[h_,\[Delta]x\:0303_]:=(t4[h,\[Delta]x\:0303]t1[h](1+\[Delta]7)+\[Delta]x\:0303)(1+\[Delta]8)*)


(* ::Input:: *)
(*t7[h_,\[Delta]x\:0303_,sk_,ck_]:=(-Hold[CorrectlyRound[sk]]t6[h,\[Delta]x\:0303](1+\[Delta]9)+t5[h,\[Delta]x\:0303,ck])(1+\[Delta]10)*)


(* ::Input:: *)
(*t8[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=(\[Delta]0 t0[h,sk,ck](1+\[Delta]11)+t7[h,\[Delta]x\:0303,sk,ck])(1+\[Delta]12)*)


(* ::Input:: *)
(*t9[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=t0[h,sk,ck](1-\[Delta]0)+t8[h,\[Delta]x\:0303,sk,ck,\[Delta]0]*)


(* ::Input:: *)
(*t10[h_,\[Delta]x\:0303_,sk_,ck_,\[Delta]0_]:=ReleaseHold[CoefficientList[Collect[t9[h,\[Delta]x\:0303,sk,ck,\[Delta]0],\[Delta]0],\[Delta]0]]*)


(* ::Input:: *)
(*t11[h_,\[Delta]x\:0303_,sk_,ck_]:=Module[{\[Delta]0,cl=t10[h,\[Delta]x\:0303,sk,ck,\[Delta]0]},cl[[1]]+uInterval cl[[2]]]*)


(* ::Input:: *)
(*cosImplementationRelativeError[x\:0303_,\[Delta]x\:0303_,xk_,sk_,ck_]:=t11[x\:0303-xk,\[Delta]x\:0303,sk,ck]/Cos[x\:0303+\[Delta]x\:0303+\[Zeta]0Interval x\:0303]-1*)


(* ::Input:: *)
(*Block[{\[Eta]=eta,\[Delta]1=d1,\[Delta]2=d2,\[Delta]3=d3,\[Delta]4=d4,\[Delta]5=d5,\[Delta]6=d6,\[Delta]7=d7,\[Delta]8=d8,\[Delta]9=d9,\[Delta]10=d10,\[Delta]11=d11,\[Delta]12=d12},CoefficientList[Collect[t9[h,\[Delta]x\:0303,sk,ck,\[Delta]0],\[Delta]0],\[Delta]0][[2]]]*)


(* ::Input:: *)
(*xi=accurateTablesXIntervals[[396]];*)


(* ::Input:: *)
(*at=accurateTables[396];*)


(* ::Input:: *)
(*Plot[{Min[cosImplementationRelativeError[x\:0303,Max[xi] \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]],Max[cosImplementationRelativeError[x\:0303,Max[xi] \[GothicU][1/2],at[[1]],at[[2]],at[[3]]]]},{x\:0303,Min[xi],Max[xi]},WorkingPrecision->30]*)


(* ::Input:: *)
(*Plot3D[{Min[cosImplementationRelativeError[x\:0303,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]],Max[cosImplementationRelativeError[x\:0303,\[Delta]x\:0303,at[[1]],at[[2]],at[[3]]]]},{x\:0303,Min[xi],Max[xi]},{\[Delta]x\:0303,-Min[xi]\[GothicU][1/2],Max[xi] \[GothicU][1/2]},RegionFunction->Function[{x\:0303,\[Delta]x\:0303},-x\:0303 \[GothicU][1/2]<\[Delta]x\:0303<x\:0303 \[GothicU][1/2]],WorkingPrecision->30,MeshShading->{{Automatic,None},{None,Automatic}},PlotStyle->{Red,Blue}]*)


(* ::Input:: *)
(*cosImplementationMaxRelativeErrorPerInterval[k_]:=*)
(*Module[*)
(*{at=accurateTables[k],xi=accurateTablesXIntervals[[k]],corners,r},*)
(*corners={{Min[xi],-Min[xi]\[GothicU][1/2]},{Min[xi],Min[xi] \[GothicU][1/2]},{Max[xi],-Max[xi]\[GothicU][1/2]},{Max[xi],Max[xi]\[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[cosImplementationRelativeError[#[[1]],#[[2]],at[[1]],at[[2]],at[[3]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{k,cosImplementationMaxRelativeErrorPerInterval[k]},{k,1,accurateTablesMaxIndex}]]*)


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
(*e=mullerE[cosImplementationMaxRelativeError,UseFMA->False];*)


(* ::Input:: *)
(*HexLiteral[CorrectlyRound[e,RoundingMode->TowardPositiveInfinity],Quotes->4]*)


(* ::Subsubsubsection::Closed:: *)
(*Proof That the Terms in h^2 \[Delta]x\:0303 and above Can Be Ignored*)


(* ::Input:: *)
(*su=Collect[Normal[Series[Sin[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*cu=Collect[Normal[Series[Cos[u],{u,0,3}]]/.u->(h+\[Delta]x),{\[Delta]x,h}]/.(\[Delta]x^n_:>0/;n>=2)*)


(* ::Input:: *)
(*Collect[TrigExpand[Cos[xk+u]]/.{Sin[u]->su,Cos[u]->cu},{\[Delta]x,h}]*)


(* ::Input:: *)
(*cosNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=(\[Delta]x\:0303 (x\:0303-xk)^2Sin[xk]/2)/Cos[x\:0303+\[Delta]x\:0303]*)


(* ::Input:: *)
(*(*cosNextTermRelativeSize[x\:0303_,\[Delta]x\:0303_,xk_]:=-(\[Delta]x\:0303 (x\:0303-xk)Cos[xk])/Cos[x\:0303+\[Delta]x\:0303]*)*)


(* ::Input:: *)
(*Plot3D[cosNextTermRelativeSize[x\:0303,\[Delta]x\:0303,at[[1]]],{x\:0303,Min[xi],Max[xi]},{\[Delta]x\:0303,-Max[xi]\[GothicU][1/2],Max[xi] \[GothicU][1/2]},RegionFunction->Function[{x\:0303,\[Delta]x\:0303},-x\:0303 \[GothicU][1/2]<\[Delta]x\:0303<x\:0303 \[GothicU][1/2]],WorkingPrecision->40,MeshShading->{{Automatic,None},{None,Automatic}},PlotRange->Full]*)


(* ::Input:: *)
(*nextTermsAtCorners[k_]:=Module[*)
(*{at=accurateTables[k],xi=accurateTablesXIntervals[[k]],corners,r},*)
(*corners={{Min[xi],-Min[xi]\[GothicU][1/2]},{Min[xi],Min[xi] \[GothicU][1/2]},{Max[xi],-Max[xi]\[GothicU][1/2]},{Max[xi],Max[xi]\[GothicU][1/2]}};*)
(*r=Block[{$MaxExtraPrecision=1000},Map[cosNextTermRelativeSize[#[[1]],#[[2]],at[[1]]]&,corners]];*)
(*Max[Abs/@r]*)
(*]*)


(* ::Input:: *)
(*ListLogPlot[Table[{k,nextTermsAtCorners[k]},{k,1,accurateTablesMaxIndex}]]*)


(* ::Input:: *)
(*cosNextTermMaxRelativeSize=Max[Table[nextTermsAtCorners[k],{k,1,accurateTablesMaxIndex}]];*)


(* ::Input:: *)
(*N[Log2[cosNextTermMaxRelativeSize]]*)


(* ::Input:: *)
(*End[]*)
