%#!/bin/rm
%university:TU Wien
 %#ok<*NOPTS>
 % close all

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  modelprops.testcase = 'pureBendingBeam'; %orderchange at lambda~.8
  %modelprops.testcase = 'cantilever';
  %modelprops.testcase = 'eccenCompressionBeam'; modelprops.ecc = 0.164669;
  %testcase = 'eccenCompressionBeam2D';
  modelprops.testcase = 'twoBeams';
<<<<<<< HEAD
  %[~,modelprops.ecc]=eccfromU(.978);
  %BpM=(modelprops.ecc)^2*0.0080678/1.319847665625e-05
  %BpGes=BpM/(1+BpM)
=======
  [~,modelprops.ecc]=eccfromU(0.5);
  BpM=(modelprops.ecc)^2*0.0080678/1.319847665625e-05
  BpGes=BpM/(1+BpM)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  
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
<<<<<<< HEAD
  %eltype = 'B31OS'; %Ti1moshenko 
=======
  %eltype = 'B31OS'; %Timoshenko 
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %eltype = 'B31OSH'; %Timoshenko 
  %eltype = 'B32' %Timoshenko 
  %eltype = 'B32H' %Timoshenko 
  eltype = 'B32OS'; %Timoshenko 
  %eltype = 'B32OSH'; %Timoshenko 
 
  
  
  % possible types of analysis
<<<<<<< HEAD
  modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  %modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=0;-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
=======
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  
  modelprops.numofelm = 2;
  
  epsil = 0.05;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
<<<<<<< HEAD
  plotfig= [14,15,-14]; %#ok<*NBRAK>
=======
  plotfig= [2,14]; %#ok<*NBRAK>
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
<<<<<<< HEAD
  modelprops.lambda = 0:epsil:20; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
=======
  modelprops.lambda = 0:epsil:.5; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %
  
  modelprops.profil.tw= 8.6e-3;
  %modelprops.forceAbaqus=true;
  modelprops.forceAbaqus=false; %default: false
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=false;
<<<<<<< HEAD
  modelprops.numofeigs=2;
=======
  modelprops.numofeigs=1;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  modelprops.allowComplex=true;
  main.closall=true;
  %main.closall=false;
  main.savefigures=true;
  %main.savefigures=false;
  %main.check=true;
  main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=true;
  
  %modelprops.sigma=-5;
  
  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
<<<<<<< HEAD
  %for x=0:9
   %ratio=x/10
   [~,ecc,zecc]=eccfromU(.1)
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,zecc);
 %end
=======
[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
modelprops=rmfield(modelprops,'forceAbaqus');
