%#!/usr/bin/env octave -q
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2023)
%author of subprograms: Johannes Kalliauer(2020-2023), Michal Malendowski (2019-2020)
%created: 2020

%% Start Trust-line-arch

%% Most important Parameters to set
% modelprops ... parameters for running Abaqus
%    modelprops.sectiondata_material_E ... defines the Youngs-Modulus
% sortType ... should the eigenvectors be sorted (currently not really working)
% forcedeig ... only plot a specific eigenvalue
% main ... structure with parameters for post-processing Abaqus-results
%   main.whichEV ... which eigenvector to use
%   main.Normierung ... how the eivenvector should be normalized


%% Recent Changes
%2023-02-16 JK: added comments for explanation
%2023-03-16 JK: removed old comments

%% Define Setting for run
 close all 
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8

  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]

  modelprops.numofelm = 20;
  
  
  epsil = 0.005;%  0.01;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  plotfig=[2,14,35];
  forcedeig = []; %1; % forced eigenvector number
  
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor =1;
  %
  
  modelprops.length=[];%19.074;% [m] 
  modelprops.sectiondata_material_E = 210e9; %[N/m^2]
  [~,modelprops.profil] =Profil('MalendowskiTLArch');
  modelprops.lambda = 0:epsil:.335;%
  modelprops.forceAbaqus=0; 
  modelprops.forcerun=1; % false... do not force it; 0.5 force if it too less lambda, 1 ... always force it.
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  main.closall=false;
  main.savefigures=1; % false... no figures, true... figures, 2 for TeX
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=false;
  main.rstabil=NaN;
  modelprops.MeterValue=1;%1000mm=1m=0.001km
  main.whichEV='bungle'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; 'Hyb'; 'bungle_rKr'; 'skip';  'sqrtK_r' 'k11' 'k0_11'
  main.Normierung='R1'; % 'skip' 'R1' 'rCT_K0_r' 'k11' 'k0_11'
  main.rho='R1'; % KtR1 R1 'skip'
  
  modelprops.followsigma=false;
  modelprops.sortJKeigval=1; %1..closest to zero, -1 ..most negative one
  main.xBezug='n'; %n..normalisiert; d..differenz zu Refwert; 1...Abaqus-Lambda; s...Stepnumber; i..individual

eltypes={'B31','B32','B32H','B33','B33H'};
for i=1:numel(eltypes)
 modelprops.elementtype = char(eltypes(i));
 main.colorshift=3*i-3;
 % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
 
end

