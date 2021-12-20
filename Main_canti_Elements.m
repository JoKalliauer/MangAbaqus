%#!
%university:TU Wien
 %#ok<*NOPTS>
 clear model
 close all
 format longG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  %modelprops.testcase = 'pureBendingBeam'; %orderchange at lambda~.8
  modelprops.testcase = 'cantilever';
  %modelprops.testcase = 'eccenCompressionBeam'; modelprops.ecc = 0.164669;
  %testcase = 'eccenCompressionBeam2D';
  
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
  eltype = 'B32OS'; %Timoshenko 
  %eltype = 'B32OSH'; %Timoshenko 
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  
  modelprops.numofelm = 20;
  
  epsil = 0.005;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  plotfig= [11,15,16]; %#ok<*NBRAK> #ok<NASGU>
  plotfig=[15,16,44,943:945,947:954];
  plotfig=[15,952,955:956];
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
  main.whichEV='bungle_rK0r'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; main.whichEV='bungle_rKr'; main.whichEV='bungle_rK0r';
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:1; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=false; %default: false
  modelprops.forcerun=0; %default=true
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  main.closall=0;
  main.savefigures=1;
  %main.savefigures=false;
  %main.check=true;
  main.check=false;
  main.colorshift=0;
  
  modelprops.sigma=-0.5;
  
%[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
modelprops=rmfield(modelprops,'forceAbaqus');

eltypes={'B31OSH','B32OSH'}%
%plotfig=[];
for i=1:numel(eltypes)
 elementtype = char(eltypes(i))
 % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 main.colorshift=i-1;
 [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,[],elementtype);
 res.stability_limit
 res.stability_limit
end
