%#!/usr/bin/env octave -q
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2023)
%author of subprograms: Johannes Kalliauer(2020-2023), Michal Malendowski (2019-2020)

%% Recent Changes
%2023-02-16 JK: added comments for explanation
%2023-05-11 JK: changed to Malendowski-dimensions
%2023-05-23 JK: close old graphs, change whichEV/Normierung, disable check

close all 
delete(findall(0,'type','figure','tag','TMWWaitbar'))
%warning('on','all')

%% Define Settings for run
modelprops.testcase = 'Kreis_2024';
modelprops.numofelm = 10;
modelprops.epsilon = 0.002;
modelprops.lambda = 0:modelprops.epsilon:.335;%.335;
modelprops.length=19.074;% [m] 
modelprops.sectiondata_material_E = 210e9; %[N/m^2]
[~,modelprops.profil] =Profil('MalendowskiTLArch');
[~,modelprops.profil] =Profil('PavlicekPage93');
modelprops.numofeigs=3;
modelprops.allowComplex=0;
modelprops.forcerun=1;
modelprops.forceAbaqus=0;
modelprops.sortJKeigval=1; %1..closest to zero, -1 ..most negative one
modelprops.elementtype = 'B32H';
main.whichEV='skip'; % 'bungle' 'NoHyb' 'skip' 'k0_11' '2023-12' '2023_12Hyb' '2023_12noHyb' '2023_12half'
main.Normierung='skip'; % 'R1' 'rNCT_K0_r' 'rCT_K0_r' 'skip 'k0_11'
main.savefigures=true;
modelprops.RefLast='old'; %'old' ; 'Wende'; 'EA'
main.xBezug='1'; %n..normalisiert; d..differenz zu Refwert; 1...Abaqus-Lambda; s...Stepnumber; i..individual
if strcmp(modelprops.RefLast,'EA')
 modelprops.epsilon=modelprops.epsilon/4000;
 modelprops.lambda=modelprops.lambda/4000;
elseif strcmp(modelprops.RefLast,'Wende')
 modelprops.epsilon=modelprops.epsilon/.15;
 modelprops.lambda=modelprops.lambda/.15;
end
main.check=0;
modelprops.alphaDRW=1;
%modelprops.alphaH=NaN;
[res,model] = Abaqus_single_run(modelprops,'none',[2,14,15,976,976.1],1,main);

