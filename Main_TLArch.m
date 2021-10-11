%#!
%university:TU Wien
 %#ok<*NOPTS>
 close all 
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  modelprops.testcase = 'TL_arch3D'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch3DKg'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch3D_sin'; %fails at ~lamdba=0.8
  %modelprops.testcase = 'TL_arch_Hinge';
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
  %eltype = 'B22'; %Timoshenko 
  %eltype = 'B22H'; %Timoshenko 
  %3D
  %eltype = 'B33'; %Euler-Bernoulli 
  %eltype = 'B33H' %Euler-Bernoulli
  %eltype = 'B31' %Timoshenko 
  %eltype = 'B31H' %Timoshenko 
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
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'Kg';
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  modelprops.numofelm = 20;
  
  
  epsil = 0.01;%  0.01;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [1,2,3,12,14,15,16,19,21,24,22,25]; %#ok<*NBRAK>
  %plotfig= [3,7,14,15,28,33]
  %plotfig= [14,28,33];
  %plotfig=[1:14,21,24,26,30,211];
  %plotfig=30;
  %plotfig=[2,7,14,21,26,211,30,34];
  %plotfig=[14,15,16,37,38,900,211];
  plotfig=[7,14,15,30,211,43]; %#ok<NASGU>
  plotfig=15;main.savefigures=1
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  maxload=30
  maxlambda=maxload/83.3;
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:maxlambda; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1;
  %
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=0; 
  modelprops.forcerun=1; % false... do not force it; 0.5 force if it too less lambda, 1 ... always force it.
  %modelprops.forcerun=false;
  modelprops.numofeigs=1;
  modelprops.allowComplex=false;
  %main.closall=true;
  main.closall=false;
  %main.savefigures=0; % false... no figures, true... figures, 2 for TeX
  %main.savefigures=false;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=false;
  main.rstabil=0.9999999960;%TL_arch3D-B31H-10-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
  %main.rstabil=0.9999999;
  modelprops.MeterValue=1;
  main.whichEV='bungle_rKr'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; main.whichEV='bungle_rKr';
  
  modelprops.followsigma=false;

  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
%[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);

%close all
% %eltypes={'B32','B32H','B31','B31H','B33','B33H'}
% % eltypes={'B32','B32H','B31','B33'} %B31H/B33H dofs aufpassen fuer rhoBungle
%eltypes={'B32','B32H','B31','B33'}
eltypes={'B32','B31','B33','B32H'}
for i=1:numel(eltypes)
 modelprops.elementtype = char(eltypes(i))
 main.colorshift=i-1;
 % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
 
end

% eltype = 'B32'
% %plotfig=902;
% numofelms = {4,8,16,32,64,128,256};
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end

% modelprops.numofelm=20;
% %% epsilon
% epsils={.2,.1,.05,.02,.01,.005,.002};
%   for i= 1:numel(epsils)
%    %modelprops.numofelm = cell2mat(i);
%    modelprops.epsilon = cell2mat(epsils(i));
%    modelprops.lambda = 0:modelprops.epsilon:max(4,20*modelprops.epsilon);
%    %    plotfig=[];
%    modelprops.ask_delete=true;
%    
%    main.colorshift=i-1;
% 
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
% 
%   end
