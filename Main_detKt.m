%#!/usr/bin/env octave
%university:TU Wien
 %#ok<*NOPTS>
 close all
 delete(findall(0,'type','figure','tag','TMWWaitbar'))

  % there are following predefined test cases:
  modelprops.testcase = 'detKt2D'; % Knickstab mitte
  %modelprops.testcase = 'd2bock';
  
  modelprops.length = 5;
  
  % possible element types (be aware of 2D and 3D):
  %2D:
  %eltype = 'B23'; %Euler-Bernoulli 
  %eltype = 'B23H'; %Euler-Bernoulli 
  eltype = 'B21'; %Timoshenko
  %eltype = 'B21H'; %Timoshenko
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
  %eltype = 'B32OSH'; %Timoshenko 
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=0;-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  
  modelprops.numofelm = 2;
  
  epsil = 1;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %sortType = 'forwardJK';
  plotfig= [900,913,914]; %#ok<*NBRAK>
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:100; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 0;
%   modelprops.MeterValue=1;
  %
  
%   modelprops.a=.15;
  
  
  %modelprops.forceAbaqus=true;
  modelprops.forceAbaqus=false; %default: false
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=false;
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  %main.closall=true;
  main.closall=false;
  %main.savefigures=true;
  main.savefigures=false;
  main.check=true;
  main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=true;
  
  %modelprops.sigma=-5;
  
  %  modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% list=[1,2,3,4,5,6,10,20,50];
list=[4,8,16,32,64,128,256];
%list=[2];
%modelprops.numelFac=[1 2];
%list=[5,50]
for i=1:numel(list)
 numofelm=list(i);
 main.colorshift=i-1;
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,numofelm);
end

modelprops=rmfield(modelprops,'forceAbaqus');
