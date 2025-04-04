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

%% Define Setting for run
 close all 
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  % there are following predefined test cases:
  modelprops.testcase = 'Kreis_arch3D';
  
  % possible element types (be aware of 2D and 3D):
  %2D:
  %eltype = 'B23'; %Euler-Bernoulli 
  %eltype = 'B22'; %Timoshenko 
  %eltype = 'B22H'; %Timoshenko 
  %3D
  %eltype = 'B33'; %Euler-Bernoulli 
  %eltype = 'B33H' %Euler-Bernoulli
  %eltype = 'B31' %Timoshenko 
  %eltype = 'B31H' %Timoshenko 
  %modelprops.elementtype = 'B32' %Timoshenko 
  %eltype = 'B32H' %Timoshenko 
  %eltype = 'B32OS'; %Timoshenko Verschiebungselement
  %eltype = 'B32OSH'; %Timoshenko hybrid
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'Kg';
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  modelprops.numofelm = 20;
  
  
  epsil = 0.005;%  0.01;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [14,28,33];
  %plotfig=[1:14,21,24,26,30,211];
  %plotfig=[2,7,14,21,26,211,30,34];
  %plotfig=[14,15,16,37,38,900,211];
  %plotfig=[14,15,16,30,43];
  %plotfig=[15,947,949,952,955:956,953];
  %plotfig=[15,45,35,19,52];%EW
  %plotfig=[45,35,19];%EW
  %plotfig=[35,19,908,916,963,919,969,902,906,917];%Verschiebungen
  %plotfig=[15,35,945,958,972];
  plotfig=[14,19,35,45,976,976.1];
  forcedeig = 1; %1; % forced eigenvector number
  
  
  %modelprops.lambda = 0:epsil:.5; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor =1;
  %
  
  modelprops.length=19.074/2;% [m] 
  modelprops.sectiondata_material_E = 210e9; %[N/m^2]
  %[~,modelprops.profil] =Profil('MalendowskiTLArch'); modelprops.lambda = 0:epsil:.335;% modelprops.numofeigs=3; % 'PavlicekPage93' 'MalendowskiTLArch'
  %[~,modelprops.profil] =Profil('PavlicekPage93'); modelprops.lambda = 0:epsil:.25; % 'PavlicekPage93' 'MalendowskiTLArch'
  modelprops.profil.h=1.54919333848297;% [m] Height of the rectangular seciton
  modelprops.profil.b=0.000161374306091976;% [m] widh of the rectangular sections
  modelprops.lambda = 0:epsil:1;
  modelprops.forceAbaqus=0; 
  modelprops.forcerun=0; % false... do not force it; 0.5 force if it too less lambda, 1 ... always force it.
  %modelprops.forcerun=false;
  modelprops.numofeigs=14;
  modelprops.allowComplex=true;
  %main.closall=true;
  main.closall=false;
  main.savefigures=1; % false... no figures, true... figures, 2 for TeX
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=false;
  %main.rstabil=0.9999999960;%TL_arch3D-B31H-10-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
  %main.rstabil=0.9999999;
  main.rstabil=NaN;
  modelprops.MeterValue=1;%1000mm=1m=0.001km
  main.whichEV='k0_11'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; 'Hyb'; 'bungle_rKr'; 'skip';  'sqrtK_r' 'k11' 'k0_11'
  main.Normierung='k0_11'; % 'skip' 'R1' 'rCT_K0_r' 'k11' 'k0_11'
  main.rho='R1'; % KtR1 R1 'skip'
  
  modelprops.followsigma=false;
  modelprops.sortJKeigval=1; %1..closest to zero, -1 ..most negative one
  main.xBezug='n'; %n..normalisiert; d..differenz zu Refwert; 1...Abaqus-Lambda; s...Stepnumber; i..individual

eltypes={'B32','B32H'};
for i=1:numel(eltypes)
 modelprops.elementtype = char(eltypes(i));
 main.colorshift=3*i-3;
 % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
 
end

% eltype = 'B32'
% %plotfig=902;
% numofelms = {4,8,16,32,64,128,256};
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end

% modelprops.numofelm=20;
% %% epsilon
% epsils={.2,.1,.05,.02,.01,.005,.002};
%   for i= 1:numel(epsils)
%    %modelprops.numofelm = cell2mat(i);
%    modelprops.epsilon = cell2mat(epsils(i));
%    modelprops.lambda = 0:modelprops.epsilon:max(4,20*modelprops.epsilon);
%    %    plotfig=[];
%    modelprops.ask_delete=true;
%    
%    main.colorshift=i-1;
% 
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
% 
%   end
