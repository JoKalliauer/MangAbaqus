<<<<<<< HEAD
%#!
%university:TU Wien
 %#ok<*NOPTS>
 % close all
 
 delete(findall(0,'type','figure','tag','TMWWaitbar'))

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch3DKg'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch3D_sin'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch_Hinge';
=======
%#!/bin/rm
%university:TU Wien
 %#ok<*NOPTS>
 % close all

  % there are following predefined test cases:
  modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  modelprops.testcase = 'TL_arch_Hinge';
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %modelprops.testcase = 'TL_arch3D_Hinge';
  %modelprops.testcase = 'pureBendingBeam'; %orderchange at lambda~.8
  %modelprops.testcase = 'cantilever';
  %modelprops.testcase = 'eccenCompressionBeam'; modelprops.ecc = 0.164669;
  %testcase = 'eccenCompressionBeam2D';
  
  %modelprops.length = [];
  modelprops.length = 5;
  
  % possible element types (be aware of 2D and 3D):
  %2D:
  %eltype = 'B23'; %Euler-Bernoulli 
<<<<<<< HEAD
  %eltype = 'B22'; %Timoshenko 
=======
  eltype = 'B22'; %Timoshenko 
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %eltype = 'B22H'; %Timoshenko 
  %3D
  %eltype = 'B33'; %Euler-Bernoulli 
  %eltype = 'B33H' %Euler-Bernoulli
  %eltype = 'B31' %Timoshenko 
  %eltype = 'B31H' %Timoshenko 
<<<<<<< HEAD
  eltype = 'B32' %Timoshenko 
=======
  %eltype = 'B31OS'; %Timoshenko 
  %eltype = 'B31OSH'; %Timoshenko 
  %eltype = 'B32' %Timoshenko 
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %eltype = 'B32H' %Timoshenko 
  %eltype = 'B32OS'; %Timoshenko 
  %eltype = 'B32OSH'; %Timoshenko 
 
  
  
  % possible types of analysis
<<<<<<< HEAD
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
=======
  modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
<<<<<<< HEAD
  %modelprops.typeofanalysis = 'Kg';
  modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  modelprops.numofelm = 20;
  numofelms = {2,5,10,20,50,100,200,500,1000,2000};
  
  epsil = 0.02;% finite difference step %epsil = 0.01;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [1,2,3,12,14,15,16,19,21,24,22,25]; %#ok<*NBRAK>
  %plotfig= [3,7,14,15,28,33]
  %plotfig= [14,28,33];
  %plotfig=[1:14,21,24,26,30,211];
  %plotfig=30;
  %plotfig=[2,7,14,21,26,211,30,34];
  plotfig=[14];
=======
  
  modelprops.numofelm = 4;
  
  epsil = 0.002;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [1,2,3,12,14,15,16,19,21,24,22,25]; %#ok<*NBRAK>
  plotfig= [14,15,211]
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:1; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %
  
  modelprops.profil.tw= 8.6e-3;
  %modelprops.forceAbaqus=true;
  modelprops.forceAbaqus=false; %default: false
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=false;
  modelprops.numofeigs=1;
<<<<<<< HEAD
  modelprops.allowComplex=false;
=======
  modelprops.allowComplex=true;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  main.closall=true;
  %main.closall=false;
  main.savefigures=true;
  %main.savefigures=false;
<<<<<<< HEAD
  main.check=true;
  %main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=false;
  main.rstabil=0.9999999960;%TL_arch3D-B31H-10-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
  %main.rstabil=0.9999999;
  
  modelprops.followsigma=false;

  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%modelprops=rmfield(modelprops,'forceAbaqus');


% for i=numofelms
%  modelprops.numofelm = cell2mat(i)
% %for j=epsils
%  %modelprops.epsilon = cell2mat(j);
%  %modelprops.lambda = 0:modelprops.epsilon:max(4)
%  plotfig=[];
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end
=======
  %main.check=true;
  main.check=false;
  main.colorshift=0;
  
  %modelprops.sigma=-5;

  % modelprops.forceAbaqus=true; modelprops.forcerun=true;
[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
modelprops=rmfield(modelprops,'forceAbaqus');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
