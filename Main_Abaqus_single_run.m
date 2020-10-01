%#!/bin/rm
%university:TU Wien
 %#ok<*NOPTS>
 % close all

  % there are following predefined test cases:
  %testcase = 'TL_arch';
  %testcase = 'TL_arch3D';
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  %testcase = 'pureBendingBeam';
  modelprops.testcase = 'cantilever';
  %testcase = 'eccenCompressionBeam';
  %     testcase = 'eccenCompressionBeam2D';
  
  %     modelprops.length = [];
  modelprops.length = 5;
  %     modelprops.ecc = [];
  %modelprops.ecc = 5;
  
  % possible element types (be aware of 2D and 3D):
  %eltype = 'B22';
  % eltype = 'B22H';
  % eltype = 'B23';
  % eltype = 'B31';
  % eltype = 'B33';
  % eltype = 'B31OS';
  % eltype = 'B31OSH';
  %eltype = 'B31'
  %eltype = 'B31H'
  %eltype = 'B32'
  %eltype = 'B32H'
  %eltype = 'B32OS';
  eltype = 'B32OSH';
  
  % possible types of analysis
  %typeofanal = 'I' %identity matrix
  %typeofanal = 'CLE' %
  %typeofanal = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3' %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4' %[ Kt0 - EW * (Kts+Ktu) ]
  
  modelprops.numofelm = 10;
  
  epsil = 0.05;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  plotfig= [1,2,12,13]; %#ok<*NBRAK>
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 5*epsil:5*epsil:10; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=true;
  %modelprops.forceAbaqus=false; %default: false
  modelprops.forcerun=true;
  %modelprops.forcerun=false; %default=true
  modelprops.closall=true;
  %modelprops.closall=false;
  %modelprops.savefigures=true;
  modelprops.savefigures=false;
  modelprops.numofeigs=11;
  modelprops.sigma=-25;
  main.check=false;
  
[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
modelprops=rmfield(modelprops,'forceAbaqus')
