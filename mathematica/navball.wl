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
smallMarking[{x,y},ToString[Abs[y/\[Pi]*180]],markings,If[y>0,n,s]]


longitude[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
ToString[Mod[x,2\[Pi]]/\[Pi]*180],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


ra[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
If[
x==0&&y==0,
"\[AriesSign]",
ToString[Mod[x,2\[Pi]]/\[Pi]*12]],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


hdg[{x_,y_},n_,s_,eq_,markings_]:=
bigMarking[
{x,y},
If[
x==0,
"N",
ToString[Mod[x,2\[Pi]]/\[Pi]*180]],
markings,
Which[y>0,n,y<0,s,y==0,eq]]


thinParallel[{x_,y_},halfLength_]:=
Style[Line[{{x-halfLength/Cos[y],y},{x+halfLength/Cos[y],y}}],Antialiasing->False]


fullParallel[y_]:=Line[{{-\[Pi],y},{\[Pi],y}}]


parallelPair[y_]:=fullParallel/@{-y,y}


parallels[markings_]:=
{markings,
Thickness[3/1024],
(*Equator*)
fullParallel[0],
(*Crosshairs along thick meridians*)
Table[
If[
Abs[y]!=\[Pi]/2&&y!=0,
thinParallel[{x,y},\[Pi]/48]],
{x,-\[Pi],\[Pi],\[Pi]/4},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}],
Table[
If[
Abs[y]<85\[Degree]&&Mod[y,\[Pi]/12]!=0,
thinParallel[{x,y},\[Pi]/120]],
{x,-\[Pi],\[Pi],\[Pi]/4},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/36}],
(*Crosshairs along thin meridians*)
Table[
If[Abs[y]<=\[Pi]/4&&Mod[x,\[Pi]/4]!=0&&Abs[y]!=\[Pi]/2&&y!=0,
thinParallel[{x,y},\[Pi]/96]],
{x,-\[Pi],\[Pi],\[Pi]/12},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}],
Table[
If[Abs[y]<=\[Pi]/4&&Mod[x,\[Pi]/4]!=0&&Abs[y]<85\[Degree]&&Mod[y,\[Pi]/12]!=0,
thinParallel[{x,y},\[Pi]/240]],
{x,-\[Pi],\[Pi],\[Pi]/12},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/36}],
parallelPair[\[Pi]/4],
parallelPair[17\[Pi]/36],
(*Polar caps*)
Rectangle[{-\[Pi],# \[Pi]/2},{\[Pi],#(\[Pi]/2-2.5\[Degree])}]&/@{-1,1}};


latitudes[n_,s_,markings_]:=
Table[
If[
y!=0&&Cos[y]!=0&&
(Mod[x,\[Pi]/4]==\[Pi]/8&&Abs[y]<75\[Degree]||
Mod[x,\[Pi]/2]==\[Pi]/4&&Abs[y]>=75\[Degree]),
latitude[{x,y},n,s,markings]],
{x,-\[Pi],\[Pi],\[Pi]/8},
{y,-\[Pi]/2,\[Pi]/2,\[Pi]/12}];


ClearAll[meridian];
meridian[\[Lambda]0_,colour_,width_,{minheight_,maxheight_}]:=
meridian[\[Lambda]0,colour,width,{minheight,maxheight}]=
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
{colour,
Style[
Rectangle[
{\[Lambda]0-ArcSin[width],minheight},
{\[Lambda]0+ArcSin[width],maxheight}],
Antialiasing->True],
Thickness[ArcSin[width]/(2\[Pi])],
Line[{{\[Lambda]0,minheight},{\[Lambda]0,maxheight}}]}]]


symmetricMeridian[\[Lambda]0_,colour_,width_,height_]:=meridian[\[Lambda]0,colour,width,{-height,height}]


fullMeridian[\[Lambda]0_,colour_,width_]:=symmetricMeridian[\[Lambda]0,colour,width,\[Pi]/2]


meridians[markings_,prime_,anti_]:=
Show[
symmetricMeridian[#,markings,.002,\[Pi]/4+\[Pi]/48]&/@Range[-11\[Pi]/12,11\[Pi]/12,\[Pi]/12],
fullMeridian[0,prime,.03],
fullMeridian[-\[Pi],anti,.03],
fullMeridian[\[Pi],anti,.03],
fullMeridian[#,markings,.01]&/@Range[-\[Pi],\[Pi],\[Pi]/4]]


background[n_,s_,eq_]:=
Graphics[
{s,Rectangle[{-\[Pi],-\[Pi]/2},{\[Pi],0}],
n,Rectangle[{-\[Pi],0},{\[Pi],\[Pi]/2}],
eq,Rectangle[{-\[Pi],-\[Pi]/36},{\[Pi],\[Pi]/36}]},
ImageMargins->0,
ImagePadding->None,
PlotRange->{{-\[Pi],\[Pi]},{-\[Pi]/2,\[Pi]/2}},
ImageSize->1024]


tightTrim=
(ImageTake[#,{1,ImageDimensions[#][[2]]-1},{2,ImageDimensions[#][[1]]}]&)@*
(ImageTrim[
#,
Last/@
Select[
Flatten[
MapIndexed[
{#1,{#2[[2]]-1,#2[[1]]-1}}&,
ImageData[#,DataReversed->True],
{2}],
1],
Function[x,Norm[x[[1]]]>.9]]]&);


(*Used in the barycentric navball*)
equatorialDisk[x_,text_,markings_]:=
{Disk[{x,0},\[Pi]/36],
Inset[
tightTrim@Rasterize[
Style[
text,
markings,
Bold,
21,
FontFamily->"Palatino"],
Background->None],
{x,0},
Center,
{Automatic,\[Pi]/36}]}


horizon=RGBColor[0.4,0.54,0.8];
sky=RGBColor[0,0.35,0.8];


navballTexture[n_,s_,markings_,eq_,prime_,anti_,longitudeFunction_,surtout_:Null]:=
Rasterize[
Show[
background[n,s,eq],
meridians[markings,prime,anti],
Graphics[
{parallels[markings],
latitudes[n,s,markings],
Table[
longitudeFunction[{x,y},n,s,eq,markings],
{x,-\[Pi],\[Pi],\[Pi]/4},
{y,-\[Pi]/4,\[Pi]/4,\[Pi]/4}],
surtout}]],
RasterSize->{1024,512}]


barycentric=navballTexture[
Lighter@Purple,
Darker@Purple,
White,
Gray,
Red,
Darker@Green,
Function[
{xy,n,s,eq,markings},
If[
Mod[xy[[1]],\[Pi]]!=0||xy[[2]]!=0,
longitude[xy,n,s,eq,markings]]],
{Purple,
equatorialDisk[0,"I",White],
equatorialDisk[-\[Pi],"II",White],
equatorialDisk[\[Pi],"II",White]}];


inertial=navballTexture[
Gray,
Black,
White,
Darker@Darker@Blue,
Red,
Darker@Green,
ra];


compass=navballTexture[
sky,
Darker@Orange,
White,
horizon,
Red,
Darker@Green,
hdg];


(* ::Code:: *)
(*SetDirectory[ParentDirectory[NotebookDirectory[]]<>"/ksp_plugin_adapter/assets"];*)
(*Export["navball_barycentric.png",barycentric];*)
(*Export["navball_inertial.png",inertial];*)
(*Export["navball_compass.png",compass];*)
