%#!/bin/rm
%university:TU Wien
 %#ok<*NOPTS>
 % close all

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  %modelprops.testcase = 'pureBendingBeam'; %orderchange at lambda~.8
  %modelprops.testcase = 'cantilever';
  %modelprops.testcase = 'eccenCompressionBeam'; modelprops.ecc = 0.164669;
  %testcase = 'eccenCompressionBeam2D';
  modelprops.testcase = 'twoBeams';
  [~,modelprops.ecc]=eccfromU(.5);
  %BpM=(modelprops.ecc)^2*0.0080678/1.319847665625e-05
  %BpGes=BpM/(1+BpM)
  
  %modelprops.length = [];
  modelprops.length = 5;
  
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
  %eltype = 'B31OS'; %Timoshenko 
  %eltype = 'B31OSH'; %Timoshenko 
  %eltype = 'B32' %Timoshenko 
  %eltype = 'B32H' %Timoshenko 
  %eltype = 'B32OS'; %Timoshenko 
  eltype = 'B32OSH'; %Timoshenko 
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=0;-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  
  modelprops.numofelm = 10;%
  
  epsil = .2;  % finite difference step %epsil = 0.005;
  sortType = 'forwardJK'; % eigenvectors sorting type: 'none', 'forwards', 'backwards',  'forwardJK'
  %plotfig= [14,15,-14,900,902,906,908]; %#ok<*NBRAK>
  plotfig=[14,15,43];
  forcedeig = [3,4]; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:20; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %
  
  modelprops.profil.tw= 8.6e-3;
  %modelprops.forceAbaqus=true;
  modelprops.forceAbaqus=0; %default: false
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=0;
  modelprops.numofeigs=4;
  modelprops.allowComplex=true;
  main.closall=true;
  %main.closall=false;
  main.savefigures=true;
  %main.savefigures=false;
  %main.check=true;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.whichEV='bungle'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; main.whichEV='bungle_rKr';
  %main.rstabil=0.9999998;
  
  %modelprops.sigma=-5;
  
  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
  %for x=0:9
   %ratio=x/10
   [~,ecc,zecc]=eccfromU(.1)
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,[],zecc);
 %end

modelprops=rmfield(modelprops,'forceAbaqus');
