%#!
%university:TU Wien
 %#ok<*NOPTS>
 %clear
 close all
% format shortG
format longG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
  %#ok<*NBRAK>
  
  
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
  %modelprops.ecc=0;
  [~,modelprops.ecc]=eccfromU(0.5);
  modelprops.testcase = 'eccenCompressionBeam'; 
  %modelprops.orientate=56;
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
  %eltypes={'B32OS','B32OSH','B31OS','B33'};
  %eltypes={'B31OS'};
  %eltypes={'B32OSH'};
  eltypes={'B32OS','B32OSH'};
 
  
  
  % possible types of analysis
  modelprops.typeofanalysis = 'I'; modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  %modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  %modelprops.numofelm = 4; %20
  
  epsil = .02; %epsil = 0.02;  % finite difference step %epsil = 0.005;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %sortType = 'forwardJK';
  %plotfig= [2,3,14,15,26,28,33]; %#ok<*NBRAK>
  %plotfig= [36,900,908,902,916,913]; %#ok<*NBRAK> 36,900,908,902,916,
  %plotfig=[0,14,36,21,211,22,18,902,2147483646,902:909,915:917] %#ok<NASGU>
  %plotfig=[7,14,15,23,30,211,36,913,908,916,906,902,908,911,916]; %#ok<NASGU>
  %plotfig=[902,211]; %#ok<NASGU>
  %plotfig=15;main.savefigures=1
  %plotfig=[902,908,916,9021,9022,913,900];
  %plotfig=[11,12,15,19,35,36,37]
  %plotfig=[957];
  plotfig=[15,45,35];
  forcedeig = []; %1; % forced eigenvector number 'none' sorting

  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:max([5.55,20*epsil]);%3.07999
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1.0;
  %int
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=0; %-1..returns error if not exist, 0..use old if exist, 1.. force new calc
  modelprops.forcerun=0; %0..use existing one, 0.5.. force run if last lambda smaller than requested, always fore a new calc.
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  main.closall=0;
  main.savefigures=1;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.rsame=0.8;
  main.rstabil=0.99999;
  main.whichEV='bungle'; % main.whichEV='bungle'; 'Disp'; 'Rot'; 'wrap'; 'Hyb'; 'bungle_rKr'; 'skip' ; 'bungle_rK0r'
  modelprops.MeterValue=100; %1000mm=1m=0.001km
  
  modelprops.sigma=0;
  modelprops.followsigma=true;
  
 %main.check=false;
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 %[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
  
% set(0, 'DefaultFigureWindowState', 'normal');

%numofelms = {1,2,4,8,16,32,64,128,256,512,1024};%numofelms = {2,5,10,20,50,100,200,500,1000,2000}
%numofelms = {3,5,6,7,9,10};
%numofelms = {1,2};
numofelms={2};
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
% end
%Exz={50,20,10,5,2,1,.5,.2,.1,.05,.02,.01,.005,0.002,0.001,.0005,0.0002,0.0001,0};
%Exz={10,1,.1,0.01,.001,.0001,0};
%Exz={0.02,0.05,.1,.2,.5,1,2,5};
%Exz={0.1,0.001,0.01,0.1,1,10,.1};
%Exz={0.05};
Exz={modelprops.ecc};
Exz={0};

%  modelprops.orientate=5.99;

for k=1:numel(Exz)
for j=1:numel(numofelms)
 modelprops.numofelm = cell2mat(numofelms(j));
 for i=1:numel(eltypes)
  %plotfig=[];
  elementtype = char(eltypes(i))
  main.colorshift=(i-1)+(j-1)*numel(eltypes);
  % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,cell2mat(Exz(k)),elementtype);
  
 end
end
end

allFigures = findall(0,'Type','figure'); % find all figures
set(allFigures,'WindowState','normal'); % set the WindowState of all figures to normal
