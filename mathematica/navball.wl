(* ::Package:: *)

ClearAll[latitude];
latitude[{x_,y_},n_,s_,markings_]:=latitude[{x,y},n,s,markings]=Inset[
hStretch[
Rasterize[
Style[" "<>ToString[Abs[y/\[Pi]*180]]<>" ",Bold,markings,72,FontFamily->"Helvetica"],
Background->If[y>0,n,s],RasterSize->100],
Round[1/Cos[y]]],
{x,y},
Center,
{Automatic,1/8}]


ClearAll[longitude];
longitude[{x_,y_},n_,s_,eq_,markings_]:=longitude[{x,y},n,s,eq,markings]=Inset[
hStretch[
Rasterize[
Style[" "<>ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*180]<>" ",Bold,markings,72,FontFamily->"Helvetica"],
Background->If[y>0,n,If[y<0,s,eq]],RasterSize->100],
Round[1/Cos[y]]],
{x,y},
Center,
{Automatic,1/6.4}]


ClearAll[ra];
ra[{x_,y_},n_,s_,eq_,markings_]:=ra[{x,y},n,s,eq,markings]=Inset[
hStretch[
Rasterize[
Style[If[x==\[Pi]&&y==0," \[AriesSign] "," "<>ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*12]<>" "],Bold,markings,72,FontFamily->"Helvetica"],
Background->If[y>0,n,If[y<0,s,eq]],RasterSize->100],
Round[1/Cos[y]]],
{x,y},
Center,
{Automatic,1/6.4}]


ClearAll[hdg];
hdg[{x_,y_},n_,s_,eq_,markings_]:=hdg[{x,y},n,s,eq,markings]=Inset[
hStretch[
Rasterize[
Style[If[x==\[Pi]," N "," "<>ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*180]<>" "],Bold,markings,72,FontFamily->"Helvetica"],
Background->If[y>0,n,If[y<0,s,eq]],RasterSize->100],
Round[1/Cos[y]]],
{x,y},
Center,
{Automatic,1/6.4}]


lines[markings_]:={markings,
Thickness[2/1024],
Table[Line[{{x,-\[Pi]/2},{x,\[Pi]/2}}],{x,\[Pi]/4,7\[Pi]/4,\[Pi]/4}],
Thickness[3/1024],
Table[If[Abs[y]!=\[Pi]/2&&y!=0,Style[Line[{{x-\[Pi]/48/Cos[y],y},{x+\[Pi]/48/Cos[y],y}}],Antialiasing->False]],{x,0,2\[Pi],\[Pi]/4},{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}],
Table[If[Abs[y]<85\[Degree]&&Mod[y,\[Pi]/12]!=0,Style[Line[{{x-\[Pi]/120/Cos[y],y},{x+\[Pi]/120/Cos[y],y}}],Antialiasing->False]],{x,0,2\[Pi],\[Pi]/4},{y,-\[Pi]/2,\[Pi]/2,\[Pi]/36}],
Line[{{0,0},{2\[Pi],0}}],
Table[If[Abs[y]<=\[Pi]/4&&Mod[x,\[Pi]/4]!=0&&Abs[y]!=\[Pi]/2&&y!=0,Style[Line[{{x-\[Pi]/96/Cos[y],y},{x+\[Pi]/96/Cos[y],y}}],Antialiasing->False]],{x,0,2\[Pi],\[Pi]/12},{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}],
Table[If[Abs[y]<=\[Pi]/4&&Mod[x,\[Pi]/4]!=0&&Abs[y]<85\[Degree]&&Mod[y,\[Pi]/12]!=0,Style[Line[{{x-\[Pi]/240/Cos[y],y},{x+\[Pi]/240/Cos[y],y}}],Antialiasing->False]],{x,0,2\[Pi],\[Pi]/12},{y,-\[Pi]/2,\[Pi]/2,\[Pi]/36}],
Line[{{0,#},{2\[Pi],#}}]&/@{-\[Pi]/4,\[Pi]/4},
Line[{{0,# 17\[Pi]/36},{2\[Pi],# 17\[Pi]/36}}]&/@{-1,1},
Rectangle[{0,# \[Pi]/2},{2\[Pi],#(\[Pi]/2-2.5\[Degree])}]&/@{-1,1}};


latitudes15[n_,s_,markings_]:=Table[
If[
y!=0&&Cos[y]!=0&&
(Mod[x,\[Pi]/4]==\[Pi]/8&&Abs[y]<75\[Degree]||
Mod[x,\[Pi]/2]==\[Pi]/4&&Abs[y]>=75\[Degree]),
latitude[{x,y},n,s,markings]],
{x,0,2\[Pi],\[Pi]/8},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}];


tightTrim=(ImageTake[#,{1,ImageDimensions[#][[2]]-1},{2,ImageDimensions[#][[1]]}]&)@*(ImageTrim[#,Last/@Select[Flatten[MapIndexed[{#1,{#2[[2]]-1,#2[[1]]-1}}&,ImageData[#,DataReversed->True],{2}],1],Function[x,Norm[x[[1]]]>.9]]]&);


backgroundAndMeridians[n_,s_,eq_,markings_,prime_,anti_]:=
{Graphics[
{s,Rectangle[{0,-\[Pi]/2},{2\[Pi],0}],
n,Rectangle[{0,0},{2\[Pi],\[Pi]/2}],
eq,Rectangle[{0,-\[Pi]/36},{2\[Pi],\[Pi]/36}],
prime,
Rectangle[{\[Pi]-ArcSin[.03],-\[Pi]/2},{\[Pi]+ArcSin[.03],\[Pi]/2}],
anti,
Rectangle[{0,-\[Pi]/2},{ArcSin[.03],\[Pi]/2}],
Rectangle[{2\[Pi]-ArcSin[.03],-\[Pi]/2},{2\[Pi],\[Pi]/2}],
markings,
Table[Rectangle[{x-ArcSin[.01],-\[Pi]/2},{x+ArcSin[.01],\[Pi]/2}],{x,0,2\[Pi],\[Pi]/4}],
Table[Line[{{x,-\[Pi]/4},{x,\[Pi]/4}}],{x,\[Pi]/12,23\[Pi]/12,\[Pi]/12}]},
ImageMargins->0,
ImagePadding->None,
PlotRange->{{0,2\[Pi]},{-\[Pi]/2,\[Pi]/2}},
ImageSize->1024],
Plot[
{-ArcCos[.005Sec[\[Lambda]+\[Pi]/4]],ArcCos[.005Sec[\[Lambda]+\[Pi]/4]],ArcCos[-.005Sec[\[Lambda]+\[Pi]/4]],-ArcCos[-.005Sec[\[Lambda]+\[Pi]/4]]},
{\[Lambda],0,2\[Pi]},
PlotRange->{-\[Pi]/2,\[Pi]/2},
PlotStyle->Directive[markings,AbsoluteThickness[.1]],
Filling->{1->Bottom,2->Top,3->Top,4->Bottom},
FillingStyle->markings,
Axes->None],
Plot[
{-ArcCos[.005Sec[\[Lambda]-\[Pi]/4]],ArcCos[.005Sec[\[Lambda]-\[Pi]/4]],ArcCos[-.005Sec[\[Lambda]-\[Pi]/4]],-ArcCos[-.005Sec[\[Lambda]-\[Pi]/4]]},
{\[Lambda],0,2\[Pi]},
PlotRange->{-\[Pi]/2,\[Pi]/2},
PlotStyle->Directive[markings,AbsoluteThickness[.1]],
Filling->{1->Bottom,2->Top,3->Top,4->Bottom},
FillingStyle->markings,
Axes->None],
Plot[
{-ArcCos[.005Csc[\[Lambda]-#]],ArcCos[.002Csc[\[Lambda]-#]],ArcCos[-.002Csc[\[Lambda]-#]],-ArcCos[-.002Csc[\[Lambda]-#]]},
{\[Lambda],0,2\[Pi]},
PlotRange->{-\[Pi]/4-\[Pi]/48,\[Pi]/4+\[Pi]/48},
PlotStyle->Directive[markings,AbsoluteThickness[.1]],
Filling->{1->-\[Pi]/4-\[Pi]/48,2->\[Pi]/4+\[Pi]/48,3->\[Pi]/4+\[Pi]/48,4->-\[Pi]/4-\[Pi]/48},
FillingStyle->markings,
Axes->None]&/@Range[\[Pi]/12,23\[Pi]/12,\[Pi]/12],
Plot[
{-ArcCos[.03Csc[\[Lambda]]],ArcCos[.03Csc[\[Lambda]]],ArcCos[-.03Csc[\[Lambda]]],-ArcCos[-.03Csc[\[Lambda]]]},
{\[Lambda],\[Pi]/2,3\[Pi]/2},
PlotRange->{-\[Pi]/2,\[Pi]/2},
PlotStyle->Directive[prime,AbsoluteThickness[.1]],
Filling->{1->Bottom,2->Top,3->Top,4->Bottom},
FillingStyle->prime,
Axes->None],
Plot[
{-ArcCos[.03Csc[\[Lambda]]],ArcCos[.03Csc[\[Lambda]]],ArcCos[-.03Csc[\[Lambda]]],-ArcCos[-.03Csc[\[Lambda]]]},
{\[Lambda],#,#+\[Pi]/2},
PlotRange->{-\[Pi]/2,\[Pi]/2},
PlotStyle->Directive[anti,AbsoluteThickness[.1]],
Filling->{1->Bottom,2->Top,3->Top,4->Bottom},
FillingStyle->anti,
Axes->None]&/@{0,3\[Pi]/2},
Plot[
{-ArcCos[.01Sec[\[Lambda]-#]],ArcCos[.01Sec[\[Lambda]-#]],ArcCos[-.01Sec[\[Lambda]-#]],-ArcCos[-.01Sec[\[Lambda]-#]]},
{\[Lambda],0,2\[Pi]},
PlotRange->{-\[Pi]/2,\[Pi]/2},
PlotStyle->Directive[markings,AbsoluteThickness[.1]],
Filling->{1->Bottom,2->Top,3->Top,4->Bottom},
FillingStyle->markings,
Axes->None]&/@Range[0,2\[Pi],\[Pi]/4]}
