(* ::Package:: *)

hStretch[img_,\[Alpha]_]:=
ImageResize[img,{\[Alpha] ImageDimensions[img][[1]],ImageDimensions[img][[2]]}]


ClearAll[markingImage];
markingImage[text_,colour_,background_,y_]:=
markingImage[text,colour,background,y]=
hStretch[
Rasterize[
Style[" "<>text<>" ",Bold,colour,72,FontFamily->"Helvetica"],
Background->background,
RasterSize->100],
Round[1/Cos[y]]]


marking[{x_,y_},text_,colour_,background_,size_]:=
Inset[
markingImage[text,colour,background,y],
{x,y},
Center,
{Automatic,size}]


bigMarking[{x_,y_},text_,colour_,background_]:=
marking[{x,y},text,colour,background,1/6.4]


smallMarking[{x_,y_},text_,colour_,background_]:=
marking[{x,y},text,colour,background,1/8]


latitude[{x_,y_},n_,s_,markings_]:=
latitude[{x,y},n,s,markings]=
smallMarking[{x,y},ToString[Abs[y/\[Pi]*180]],markings,If[y>0,n,s]]


longitude[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*180],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


ra[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
If[
x==\[Pi]&&y==0,
"\[AriesSign]",
ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*12]],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


hdg[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
If[
x==\[Pi],
"N",
ToString[Mod[(x-\[Pi]),2\[Pi]]/\[Pi]*180]],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


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


meridian[\[Lambda]0_,colour_,width_,{minheight_,maxheight_}]:=
Show[
Plot[
{-ArcCos[width Csc[\[Lambda]-\[Lambda]0]],ArcCos[width Csc[\[Lambda]-\[Lambda]0]],ArcCos[-width Csc[\[Lambda]-\[Lambda]0]],-ArcCos[-width Csc[\[Lambda]-\[Lambda]0]]},
{\[Lambda],\[Lambda]0-\[Pi]/2,\[Lambda]0+\[Pi]/2},
PlotRange->{minheight,maxheight},
PlotStyle->Directive[colour,AbsoluteThickness[0]],(*needed for antialiasing*)
Filling->{1->minheight,2->maxheight,3->maxheight,4->minheight},
FillingStyle->colour,
Axes->False],
Graphics[
Rectangle[
{\[Lambda]0-ArcCos[width],minheight},
{\[Pi]+ArcSin[width],maxheight}]]]


fullMeridian[\[Lambda]0_,colour_,width_]:=meridian[\[Lambda]0,colour,width,{-\[Pi]/2,\[Pi]/2}]


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
meridian[#,markings,.002,{-\[Pi]/4-\[Pi]/48,\[Pi]/4+\[Pi]/48}]&/@Range[\[Pi]/12,23\[Pi]/12,\[Pi]/12],
fullMeridian[\[Pi],prime,.03],
fullMeridian[0,anti,.03],
fullMeridian[2\[Pi],anti,.03],
fullMeridian[#,markings,.01]&/@Range[0,2\[Pi],\[Pi]/4]}
