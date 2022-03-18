(* ::Package:: *)

(* ::Input:: *)
(*Quit[];*)


(* ::Subsubsection:: *)
(*Import Theory DD cs calculation and Darwin and X1T experimental constraint and Pico SD constraint*)


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
<<"Import/dd_module.wl";


(* ::Section::Closed:: *)
(*Define Functions that determine the constraint on couplings gDM from SI DD - These are rudimentary, please modify and use wisely*)


(* ::Input::Initialization:: *)
(*Tolerance, can't go all the way down to zero dark matter mass*)
mchilim=10^-2;
(******************************)
(******************************)
(*XENON1T CONSTRAINTS*)
(******************************)
(*For the uR model*)
fr3RGESIuR[msq_,mchi_]:=Module[{ggdm={1,0,0,1,0,1},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>mchilim,tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],{x,0.1,1*10^-10,100},AccuracyGoal->5],tt={x-> 100.}];x/.tt];
(*For the dR model*)
fr3RGESIdR[msq_,mchi_]:=Module[{ggdm={0,1,1,0,1,0},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],{x,0.1,1*10^-10,100},AccuracyGoal->5],tt={x-> 100.}];x/.tt];
(*For the qL model*)
fr3RGESIql[msq_,mchi_]:=Module[{ggdm={1,1,1,1,1,1},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],{x,0.1,1*10^-10,100},AccuracyGoal->5],tt={x-> 100.}];x/.tt];
(******************************)
(*XENON1T CONSTRAINTS*)
(******************************)
(*These work better in the compressed region*)
(*For the uR model*)
fr3RGESIdRv2[msq_,mchi_]:=Module[{ggdm={0,1,1,0,1,0},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];
(*For the dR model*)
fr3RGESIuRv2[msq_,mchi_]:=Module[{ggdm={1,0,0,1,0,1},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];
(*For the qL model*)
fr3RGESIqlv2[msq_,mchi_]:=Module[{ggdm={1,1,1,1,1,1},tt,r1},
If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[x1SIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];

(******************************)
(*Repeat the same for DARWIN CONSTRAINTS*)
(******************************)
(*These work better in the compressed region*)
(*For the uR model*)
fr3RGESIdRv2D[msq_,mchi_]:=Module[{ggdm={0,1,1,0,1,0},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[darwinSIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];
(*For the dR model*)
fr3RGESIuRv2D[msq_,mchi_]:=Module[{ggdm={1,0,0,1,0,1},tt},If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[darwinSIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];
(*For the qL model*)
fr3RGESIqlv2D[msq_,mchi_]:=Module[{ggdm={1,1,1,1,1,1},tt,r1},
If[mchi<mchilim,tt={x-> 100.}];If[msq>mchi  &&mchi>= mchilim,If[Abs[(1-mchi/msq)]>0.8,r1={x,0.1,0,100},r1={x,1*10^-10,1*10^-10,4 Pi}];tt=FindRoot[darwinSIdtExtrap2[mchi]==Log[sigSIRGE[msq,mmL,mchi,\[Alpha]srun[2],x*ggdm]],r1],tt={x-> 100.}];x/.tt];




(* ::Input:: *)
(*(* msq : Mediator mass, mchi: Dark matter mass.*)*)
(*msqTest=500;*)
(*mchiTest=200*)
(*(*Returns the upper bound on the coupling gDM for the uR model*)*)
(*fr3RGESIuRv2D[msqTest,mchiTest]*)
(*(*Returns the upper bound on the coupling gDM for the uR model*)*)
(*fr3RGESIdRv2D[msqTest,mchiTest]*)
(*(*Returns the upper bound on the coupling gDM for the qL model*)*)
(*fr3RGESIqlv2D[msqTest,mchiTest]*)


(* ::Section:: *)
(*Pico-60 Spin Dependent - Coming Soon*)


(* ::Section:: *)
(*Darwin Spin Independent  - Coming Soon*)


(* ::Section:: *)
(*Darwin Spin Dependent - Coming Soon*)
