(* ::Package:: *)

(* ::Section:: *)
(*Utility Functions*)


(* ::Input:: *)
(*CheckFile[filename_String]:=Module[*)
(*{oldNames,newNames,commonNames,oldValues,newValues,sameNames, sameValues},*)
(*If[Context[]!="Global`",End[]];(*Clear any leftover context.*)*)
(*Remove["old`*"];*)
(*Remove["new`*"];*)
(*Begin["old`"];*)
(*Get[filename<>".wl"];*)
(*End[];*)
(*Print[Names["old`*"]];*)
(*Begin["new`"];*)
(*Get[filename<>"_new.wl"];*)
(*End[];*)
(*Print[Names["new`*"]];*)
(*oldNames=StringReplace[Names["old`*"],"old`"->""];*)
(*newNames=StringReplace[Names["new`*"],"new`"->""];*)
(*sameNames=oldNames==newNames;*)
(*commonNames=Intersection[oldNames,newNames];*)
(*oldValues=Map[Symbol,Names[Map["old`"<>#&,commonNames]]];*)
(*newValues=Map[Symbol,Names[Map["new`"<>#&,commonNames]]];*)
(*sameValues=MapThread[#1==#2&,{oldValues,newValues}];*)
(*{sameNames,AllTrue[sameValues,TrueQ]}*)
(*]*)


(* ::Section:: *)
(*Directory Paths*)


(* ::Input:: *)
(*directory="C:\\Users\\phl\\Projects\\GitHub\\Principia\\Principia\\"*)


(* ::Input:: *)
(*mathematicaDirectory=directory<>"mathematica\\"*)


(* ::Input:: *)
(*tempDirectory=directory<>"temp\\Release\\x64\\"*)


(* ::Section:: *)
(*Astronomy Tests*)


(* ::Text:: *)
(*Run tests with:*)


 .\Release\x64\astronomy_tests.exe --gtest_also _run _disabled _tests --gtest_filter=-*SECULAR*


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"lunar_orbit_10x10.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"lunar_orbit_20x20.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"lunar_orbit_25x25.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"lunar_orbit_30x30.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"ksp_system_convergence.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"fully_perturbed_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"j2_perturbed_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"unperturbed_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"L74_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"L94_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"E01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"E14_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"G01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"L01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"R01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"C01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"C06_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"C34_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"J01_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"J07_elements.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"solar_system_convergence.generated"]*)


(* ::Input:: *)
(*CheckFile[mathematicaDirectory<>"\:043c\:043e\:043b\:043d\:0438\:044f_orbit.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"ksp_system.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"corrected.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"stock.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"phobos.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"trappist_alignments.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"trappist_transits.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"trappist_periods.generated"]*)


(* ::Section:: *)
(*Integrator Tests*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"convergence.Quinlan1999Order8A.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"convergence.Quinlan1999Order8B.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"convergence.QuinlanTremaine1990Order8.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"convergence.QuinlanTremaine1990Order10.generated"]*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"convergence.QuinlanTremaine1990Order12.generated"]*)


(* ::Section:: *)
(*Physics Tests*)


(* ::Input:: *)
(*CheckFile[tempDirectory<>"discrete_trajectory_compression.generated"]*)



