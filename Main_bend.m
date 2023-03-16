%#!/usr/bin/env octave -q
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2023)
%author of subprograms: Johannes Kalliauer(2020-2023), Michal Malendowski (2019-2020)
%created: 2020

%% Start pure-bending-beam on two supports

%% Most important Parameters to set
% modelprops ... parameters for running Abaqus
% sortType ... should the eigenvectors be sorted (currently not really working)
% forcedeig ... only plot a specific eigenvalue
% main ... structure with parameters for post-processing Abaqus-results
%   main.whichEV ... which eigenvector to use
%   main.Normierung ... how the eivenvector should be normalized


%% Recent Changes
%2023-02-21 JK: added explanatation
%2023-03-16 JK: removed old comments

%% Define Setting for run


 %#ok<*NOPTS>
 clear model
 close all
 format longG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  modelprops.testcase = 'pureBendingBeamJK'; %orderchange at lambda~.8
  
  modelprops.length = 5;%m
  
  eltypes={'B31','B31H','B32','B32H','B32OS','B32OSH','B33','B33H'};%  eltypes={'B31','B31H','B31OS', 'B31OSH','B32','B32H','B32OS', 'B32OSH','B33','B33H'};
 
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  
  modelprops.numofelm = 20;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  plotfig=[14,35];
  
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
  
  modelprops.loadfactor = 1;
  %
  
  modelprops.forceAbaqus=false; %-1 ... don't allow reruning, false... dont force rerun, 0.5 rerun if too less lambda, 1 force rerun
  modelprops.forcerun=true; %0 dont force, 0.5 force run if last lambda smaller than requested; 1 force run
  modelprops.numofeigs=1;
  modelprops.allowComplex=false;
  main.closall=0;
  main.savefigures=false; % false.. dont safe figures(faster), true safe figures (slow)
  main.check=false;
  modelprops.ask_delete=true;
  modelprops.MeterValue=1; %1000mm=1m=0.001km
  main.whichEV='bungle'; % main.whichEV='bungle'; 'Disp'; 'Rot'; 'wrap'; 'Hyb'; 'rNCT_K0_r';'rCT_K0_r'; 'split'; 'corrected' ; 'k11';  'sqrtK_r'; 'sqrtK0_r'; 'NoHyb' 'k0_11'
  main.Normierung='R1'; % 'R1'; 'rCT_K0_r'; 'A0R1'; 'sqrtK_r' 'k0_11'
  main.rho='R1'; % KtR1 R1; 'A0R1'
  main.xBezug='n'; %n..normalisiert; d..differenz zut Refwert
   
  modelprops.followsigma=true;
  
  
epsils={.05,.02,.01}
for j=1:numel(epsils)
    modelprops.epsilon = cell2mat(epsils(j));
    modelprops.lambda = 0:modelprops.epsilon:max(.8,30*modelprops.epsilon);
 for i=1:numel(eltypes)
  elementtype = char(eltypes(i));
  main.colorshift=(i-1)+(j-1)*numel(eltypes);
  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,[],elementtype);
  res.stability_limit
 end
end
  

% close all
cfig = containers.Map;
cfig('outputformats')='sep'%svg,eps,pdf

numofelm=modelprops.numofelm;
 xdata=(0:numofelm)*modelprops.length/numofelm
 
for Dir=1:6
 ydata=model.eigvecDR{20}(Dir,1:numofelm+1,3)
 myylabel=strcat('Eigenvectorcomponent',num2str(Dir));
 plotitJK(xdata,ydata,'./','x-axis of beam [m]',myylabel,myylabel,cfig,Dir+100) 
end