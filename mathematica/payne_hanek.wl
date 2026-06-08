(* ::Package:: *)

(* ::Input:: *)
(*Get[FileNameJoin[NotebookDirectory[],"ieee754_floating_point.wl"]]*)


(* ::Input:: *)
(*On[Assert]*)


(* ::Input:: *)
(*\[Alpha]=4/\[Pi]*)


(* ::Input:: *)
(*BaseForm[N[\[Alpha],100],2]*)


(* ::Text:: *)
(*Bits start at 0 and are counted with negative indices.*)


(* ::Input:: *)
(*chunk[x_,length_,start_]:=*)
(*Module[{digits=RealDigits[x,2,length,start]},*)
(*{FromDigits[digits[[1]],2],digits[[2]]-length}] *)


(* ::Input:: *)
(*Assert[chunk[\[Alpha],7,0]=={2^^1010001,-6}]*)


(* ::Input:: *)
(*Assert[chunk[\[Alpha],7,-7]=={2^^0111110,-13}]*)


(* ::Input:: *)
(*Assert[chunk[\[Alpha],7,-14]=={2^^0110000,-20}]*)


(* ::Input:: *)
(*Assert[chunk[\[Alpha],7,-21]=={2^^0110110,-27}]*)


(* ::Input:: *)
(*Assert[chunk[\[Alpha],7,-28]=={2^^1110010,-34}]*)


(* ::Input:: *)
(*Log[2,4.94``400*^-324]*)


(* ::Input:: *)
(*bitsPerChunk=Floor[binary64[[1]]/2]*)


(* ::Input:: *)
(*smallestExponent=2^(binary64[[2]]-1)+binary64[[1]]-3*)


(* ::Input:: *)
(*chunks=Table[chunk[\[Alpha],bitsPerChunk,-bitsPerChunk i],{i,0,Floor[smallestExponent/bitsPerChunk]}]*)


(* ::Input:: *)
(*hexFloatLiteral[x_Integer,exponent_Integer]:=*)
(* Module[*)
(*  {group=4,*)
(*digits=7,*)
(*groups},*)
(*groups=Ceiling[digits/group];*)
(*  StringJoin["0x",StringDrop[StringRiffle[StringPartition[IntegerString[x,16,group*groups],UpTo[group]],"'"],group*groups-digits],".0p",*)
(*If[exponent<0,"-",""],IntegerString[exponent]]]*)


(* ::Input:: *)
(*hexFloatLiteral[42722829,-25]*)


(* ::Input:: *)
(*N[16^^28be60d.0 2^-25]*)


(* ::Text:: *)
(*Note that the last chunk is not exactly representable because of progressive underflow.*)


(* ::Input:: *)
(*chunkStrings=MapIndexed[hexFloatLiteral[#1 [[1]],#1[[2]]]&,chunks]*)


(* ::Input:: *)
(*Export[*)
(*FileNameJoin[{NotebookDirectory[],"..\\numerics\\payne_hanek.mathematica.h"}],*)
(*"#pragma once*)
(**)
(*#include <array>*)
(*#include <cstdint>*)
(**)
(*namespace principia {*)
(*namespace numerics {*)
(*namespace _payne_hanek {*)
(*namespace internal {*)
(**)
(*constexpr std::int64_t PayneHanekBitsPerChunk = "<>ToString[bitsPerChunk]<>";*)
(**)
(*// Chunks of "<>ToString[bitsPerChunk]<>" bits of 4/\[Pi], up to exponent "<>ToString[-smallestExponent]<>".*)
(*constexpr std::array<double, "<>ToString[Length[chunkStrings]]<>"> PayneHanekChunks{\n    "<>*)
(*StringRiffle[chunkStrings,",\n    "]*)
(*<>"};*)
(**)
(*}  // namespace internal*)
(**)
(*using internal::PayneHanekBitsPerChunk;*)
(*using internal::PayneHanekChunks;*)
(**)
(*}  // namespace _payne_hanek*)
(*}  // namespace numerics*)
(*}  // namespace principia*)
(*",*)
(*"text"]*)
