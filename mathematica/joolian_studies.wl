(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
joolSystem={{-3*^8,3*^8},{-3*^8,3*^8},{-7*^7,7*^7}};
orthogonalOrbitPlot[q_,range_]:=TableForm[
{{ListPlot[#[[;;,{1,3}]]&/@q\[Transpose],AspectRatio->Automatic,ImageSize->{600,240},PlotRange->range[[{1,3}]],PlotStyle->Thickness[0],Joined->True],
 Graphics[
  {Line[{{-8,1},{-4,2},{-4,-2},{-8,-1},{-8,1}}],
   Circle[],Circle[{0,0},2],
   Dashing[{Medium,Tiny,Tiny,Tiny}],
   Line[{{-9,0},{3,0}}],
   Line[{{0,3},{0,-3}}]}]},
 {ListPlot[#[[;;,;;2]]&/@q\[Transpose],AspectRatio->Automatic,ImageSize->{600,600},PlotRange->range[[;;2]],PlotStyle->Thickness[0],Joined->True],
   ListPlot[#[[;;,{3,2}]]&/@q\[Transpose],AspectRatio->Automatic,ImageSize->{240,600},PlotRange->range[[{3,2}]],PlotStyle->Thickness[0],Joined->True]}}];


<<"retrobop_predictable_years.generated.wl";


Export["retrobop_1_a.png",orthogonalOrbitPlot[Reverse/@barycentricPositions1,joolSystem],ImageResolution->300];


Export["retrobop_2_a.png",orthogonalOrbitPlot[Reverse/@barycentricPositions2,joolSystem],ImageResolution->300];


Export["retrobop_5_a.png",orthogonalOrbitPlot[Reverse/@barycentricPositions5,joolSystem],ImageResolution->300];


<<"retrobop_century.generated.wl";


Export[
"retrobop_apsides.png",
ListPlot[
{#[[1]]/(60*60*24*365.25),#[[2]]}&/@#&/@{
 {polTimes,polSeparations}\[Transpose],
 {bopTimes,bopSeparations}\[Transpose],
 {tyloTimes,tyloSeparations}\[Transpose],
 {vallTimes,vallSeparations}\[Transpose],
 {laytheTimes,laytheSeparations}\[Transpose]},
PlotRange->{Full,{0,2.2*^8}},ImageSize->800,AxesLabel->{"t (a)","apsides (m)"}]]


Export[
"retrobop_eccentricities.png",
ListPlot[
{Range[Length[#]]/(24*365.25),#}\[Transpose]&/@{bopEccentricities,bopJacobiEccentricities},
PlotRange->Full,ImageSize->800,AxesLabel->{"t (a)","\!\(\*SubscriptBox[\(e\), \(Bop\)]\)"}]]


Export[
"retrobop_nodes.png",
ListPlot[
{Range[Length[#]]/(24*365.25),#}\[Transpose]&/@{bopNodes,bopJacobiNodes},
PlotRange->Full,ImageSize->800,AxesLabel->{"t (h)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(Bop\)]\) (\[Degree])"}]]


Export[
"retrobop_arguments.png",
ListPlot[
{Range[Length[#]]/(24*365.25),#}\[Transpose]&/@{bopArguments,bopJacobiArguments},
PlotRange->Full,ImageSize->800,AxesLabel->{"t (h)","\!\(\*SubscriptBox[\(\[Omega]\), \(Bop\)]\) (\[Degree])"}]]


Export[
"retrobop_apsides_1_a.png",
ListPlot[
{#[[1]]/(60*60*24*365.25),#[[2]]}&/@#&/@{
 {polTimes,polSeparations}\[Transpose],
 {bopTimes,bopSeparations}\[Transpose],
 {tyloTimes,tyloSeparations}\[Transpose],
 {vallTimes,vallSeparations}\[Transpose],
 {laytheTimes,laytheSeparations}\[Transpose]},
PlotRange->{{0,1},{0,2.2*^8}},ImageSize->800,AxesLabel->{"t (a)","apsides (m)"}]]


Export[
"retrobop_apsides_10_a.png",
ListPlot[
{#[[1]]/(60*60*24*365.25),#[[2]]}&/@#&/@{
 {polTimes,polSeparations}\[Transpose],
 {bopTimes,bopSeparations}\[Transpose],
 {tyloTimes,tyloSeparations}\[Transpose],
 {vallTimes,vallSeparations}\[Transpose],
 {laytheTimes,laytheSeparations}\[Transpose]},
PlotRange->{{0,10},{0,2.2*^8}},ImageSize->800,AxesLabel->{"t (a)","apsides (m)"}]]
