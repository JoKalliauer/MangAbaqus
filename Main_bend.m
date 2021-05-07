%#!
%university:TU Wien
 %#ok<*NOPTS>
 close all
%format shortG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  %modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %testcase = 'TL_arch_Hinge';
  %testcase = 'TL_arch3D_Hinge';
  modelprops.testcase = 'pureBendingBeamJK'; %orderchange at lambda~.8
  %modelprops.testcase = 'pureBendingBeamMalendowski';
  %modelprops.testcase = 'cantilever';
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
  eltypes={'B32OS','B32OSH','B31OS','B31OSH'}
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';%modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';%modelprops.sigma=pi() %
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
  numofelms = {2,5,10,20,50,100,200,500,1000,2000};
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [2,7,14,15,21,211,26,28,29]; %#ok<*NBRAK>
  plotfig=[14,15,16,37];
  
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  modelprops.epsilon = .02;  % finite difference step %epsil = 0.005;
  epsils= {1,.5,.2,.1,.05,.02,.01,.005,.002,.001};
  modelprops.lambda = 0:modelprops.epsilon:max(.7,2*modelprops.epsilon);%10; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.loadfactor = 1;
  %
  
  %modelprops.profil.tw= 8.6e-3;
  %modelprops.forceAbaqus=true;
  %%%%%modelprops.forceAbaqus=0.5; % force run if last lambda smaller than requested
  modelprops.forceAbaqus=false; %default: false
  %modelprops.forceAbaqus=-1; %don't run simulation even if not existing
  %modelprops.forcerun=true; %default=true
  modelprops.forcerun=0.5; % force run if last lambda smaller than requested
  %modelprops.forcerun=false;
  modelprops.numofeigs=1;
  modelprops.allowComplex=false;
  %main.closall=true;
  main.closall=false;
  main.savefigures=true;
  %main.savefigures=false;
  %main.check=true;
  main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=true;
  modelprops.MeterValue=1000;
  main.whichEV='bungle'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb';
  
  %modelprops.sigma=-10;
  modelprops.followsigma=false;
  
  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);

% eltypes={'B32OS','B32OSH','B31OS','B31OSH'}
eltypes={'B32OS','B32OSH','B31OS','B31OSH'}
%plotfig=[];
for i=1:numel(eltypes)
 elementtype = char(eltypes(i))
%  if strcmp(elementtype,'B32OSH') ||  strcmp(elementtype,'B31OSH')
%   main.check=true;
%  else
%   main.check=false;
%  end
 % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 main.colorshift=i-1;
 [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,[],elementtype);
 
end
  
%   for i= {.05,.02,.01}
%    %modelprops.numofelm = cell2mat(i)
%    modelprops.epsilon = cell2mat(i);
%    modelprops.lambda = 0:modelprops.epsilon:max(4,20*modelprops.epsilon)
%    %    plotfig=[];
%    modelprops.ask_delete=true;
% 
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
% 
%   end
