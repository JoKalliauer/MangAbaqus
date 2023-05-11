%#!/usr/bin/env octave -q
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2023)
%author of subprograms: Johannes Kalliauer(2020-2023), Michal Malendowski (2019-2020)

%% Recent Changes
%2023-02-16 JK: added comments for explanation
%2023-05-11 JK: changed to Malendowski-dimensions

%% Define Settings for run
modelprops.testcase = 'Kreis_arch3D';
modelprops.numofelm = 100;
modelprops.epsilon = 0.0005;
modelprops.length=19.074/2;% [m] 
modelprops.sectiondata_material_E = 210e9; %[N/m^2]
[~,modelprops.profil] =Profil('MalendowskiTLArch');
modelprops.numofeigs=3;
[res,model] = Abaqus_single_run(modelprops,'none',35,[]);

