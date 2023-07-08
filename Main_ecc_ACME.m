%#!
%university:TU Wien
 %#ok<*NOPTS>
 %clear
 close all
% format shortG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 %set(0, 'DefaultFigureWindowStyle', 'docked');
  %#ok<*NBRAK>
   %#ok<*NASGU>
  
  
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
  %modelprops.ecc=50; %.5
  %modelprops.ecc=.005; %.5
  %modelprops.ecc=.5;
  [~,modelprops.ecc]=eccfromU(0.5);
  modelprops.testcase = 'eccenCompressionBeam'; 
  %modelprops.testcase = 'eccenCompressionBeam64';
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
  %modelprops.elementtype = 'B32OS'; %Timoshenko 
  %modelprops.elementtype = 'B32OSH'; %Timoshenko 
  %eltypes={'B33','B33H','B31','B31H','B31OS','B31OSH','B32','B32H','B32OS','B32OSH'};
  %eltypes={'B32OS','B32OSH','B31OSH','B33','B32','B31OS'}; 
  %eltypes={'B31'}%,'B32','B32OS','B31OS','B33'};
  eltypes={'B32OS','B32OSH'}
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I'; modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=-1; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  %modelprops.numofelm = 4; %20
  
  epsil = 0.1;%.01; %epsil = 0.02; 0.005  % finite difference step %epsil = 0.005;
  sortType = 'forwardJK'; % eigenvectors sorting type: 'none', 'forwards', 'backwards', 'forwardJK'
  %plotfig= [2,3,14,15,26,28,33]; %#ok<*NBRAK>
  %plotfig= [36,900,908,902,916,913]; %#ok<*NBRAK> 36,900,908,902,916,
  %plotfig=[0,14,36,21,211,22,18,902,2147483646,902:909,915:917] %#ok<NASGU>
  %plotfig=[7,14,15,23,30,211,36,913,908,916,906,902]; %#ok<NASGU>
  %plotfigProb=[947,949,944:945,948]
  %plotfig=[11,14,15,19,43,952,955:956];  
  %plotfig=[11,12,15,19,35,36,37,45];
  %plotfig=[15,19,45];
  plotfig=[14,15,43,945];
  %plotfig=[15,19,43,45];
  %plotfig=[15,16,943,953];
  %plotfig=[902,908,916,9021,9022,913,900];
  %plotfig=[15,35,945,958,972];
  forcedeig = []; %1; % forced eigenvector number 'none' sorting

  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0*epsil:epsil:max([3.4749,20*epsil]);%3.07999
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %int
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=0; %-1..returns error if not exist, 0..use old if exist, 1.. force new calc
  modelprops.forcerun=1; %0..use existing one, 0.5.. force run if last lambda smaller than requested, always fore a new calc.
  modelprops.numofeigs=4;
  modelprops.allowComplex=1;
  %main.closall=true;
  main.closall=false;
  main.savefigures=1;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.rsame=0.8;
  main.rstabil=0.99999;
  main.whichEV='NoHyb'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; main.whichEV='bungle_rKr';
  main.Normierung='R1'; % 'R1'; 'rCT_K0_r'
  main.rho='R1'; % KtR1 R1 'A0R1' 
  modelprops.MeterValue=1; %1000mm=1m=0.001km
  
  
 %main.check=false;
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 %[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
  
% set(0, 'DefaultFigureWindowState', 'normal');

%numofelms = {1,2,4,8,16,20,32,64,128,256,512,1024};%numofelms = {2,5,10,20,50,100,200,500,1000,2000}
%numofelms = {3,5,6,7,9,10};
numofelms = {4};
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end

for j=1:numel(numofelms)
 modelprops.numofelm = cell2mat(numofelms(j));
 for i=1:numel(eltypes)
  %plotfig=[];
  elementtype = char(eltypes(i))
  main.colorshift=(i-1)+(j-1)*numel(eltypes);
  % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,modelprops.ecc,elementtype);
  
 end
end
  
%     list=[1,.1,.01,10,0.001,1000,100]
%   for i=1:numel(list)
%    modelprops.MeterValue=list(i);
%      [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%   end
  
%    
%   %eltypes={'B33','B31','B31OS','B32','B32OS'}
%   list=[2,5,10,20]
% % for i=1:numel(list)
% %  numofelm=list(i);
% %    [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,numofelm);
% % end
% 
% %% epsilon
% epsils={1,.5,.2,.1,.05,.02,.01,.005,.002}%,.01,.005,.002
%   for i=1:numel(epsils)
%    %modelprops.numofelm = cell2mat(i)
%    modelprops.epsilon = cell2mat(epsils(i));
%    modelprops.lambda = 0:modelprops.epsilon:max(5);
%    main.colorshift=i-1;
%    %    plotfig=[];
% 
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
% 
%   end
% modelprops.epsilon = 0.02;
% modelprops.lambda = 0:modelprops.epsilon:max(5);
% 
% 
% numofelms = {4,8,16,32,64,128,256};%numofelms = {2,5,10,20,50,100,200,500,1000,2000}
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end

allFigures = findall(0,'Type','figure'); % find all figures
set(allFigures,'WindowState','normal'); % set the WindowState of all figures to normal
