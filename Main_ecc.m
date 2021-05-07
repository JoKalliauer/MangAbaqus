%#!/bin/rm
%university:TU Wien
 %#ok<*NOPTS>
 % close all
 format shortG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 
  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  %modelprops.testcase = 'pureBendingBeam'; %orderchange at lambda~.8
  %modelprops.testcase = 'cantilever';
  %modelprops.ecc = 0.164669;
  %modelprops.ecc=(81*sqrt(64373/403390))/800;
  %modelprops.ecc=0;
  %modelprops.ecc=0.02;
  %modelprops.ecc=.005;
  [~,modelprops.ecc]=eccfromU(0.5);
  modelprops.testcase = 'eccenCompressionBeam'; 
  %testcase = 'eccenCompressionBeam2D';
  BpM=(modelprops.ecc)^2*0.0080678/1.319847665625e-05;
  BpGes=BpM/(1+BpM);

  
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
  modelprops.elementtype = 'B32OS'; %Timoshenko 
  %eltype = 'B32OSH'; %Timoshenko 
  %eltypes={'B33','B33H','B31','B31H','B31OS','B31OSH','B32','B32H','B32OS','B32OSH'};
  eltypes={'B33','B31','B31OS','B32','B32OS'};
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I'; modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  modelprops.numofelm = 10;
  
  epsil = .02; %epsil = 0.02;  % finite difference step %epsil = 0.005;
  %sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  sortType = 'forwardJK';
  %plotfig= [2,3,14,15,26,28,33]; %#ok<*NBRAK>
  plotfig= [36,900,908,909]; %#ok<*NBRAK>
  %plotfig=[14,15,28]
  forcedeig = []; %1; % forced eigenvector number 'none' sorting

  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:max([5,20*epsil]);%3.07999
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %int
  
  modelprops.profil.tw= 8.6e-3;
  %modelprops.forceAbaqus=true; modelprops.forcerun=true;
  %modelprops.forceAbaqus=false; %default: false
  modelprops.forceAbaqus=int8(-1); % throw an error if calculation does not exist
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=false;
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  main.closall=true;
  %main.closall=false;
  main.savefigures=true;
  %main.savefigures=false;
  %main.check=true;
  main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.rsame=0.8;
  main.rstabil=0.99999;
  
  modelprops.sigma=0;
  
 %main.check=false;
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
  

%   for i=1:numel(eltypes)
%   %plotfig=[];
%    elementtype = char(eltypes(i))
%   % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,modelprops.ecc,elementtype);
% 
%   end
  
%     modelprops.elementtype = 'B31OS'; 
%   %eltypes={'B33','B31','B31OS','B32','B32OS'}
%   list=[2,5,10,20]
% for i=1:numel(list)
%  numofelm=list(i);
%    [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,numofelm);
% end
