%#!
%university:TU Wien
 %#ok<*NOPTS>
 %clear
 clear model
 close all
 format longG
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');
 %  set(0, 'DefaultFigureWindowState', 'minimized');
 %set(0, 'DefaultFigureWindowStyle', 'docked');

  % there are following predefined test cases:
  modelprops.testcase = 'pureBendingBeamJK'; %orderchange at lambda~.8
  %modelprops.testcase = 'pureBendingBeamMalendowski';
  %modelprops.testcase = 'pureBendingCantilever'; modelprops.orientate=5;
  
  %modelprops.length = [];
  modelprops.length = 5;
  
  % possible element types (be aware of 2D and 3D):
  %2D:
  %eltype = 'B23'; %Euler-Bernoulli 
  %eltype = 'B22'; %Timoshenko 
  %eltype = 'B22H'; %Timoshenko 
  %3D
  %modelprops.elementtype = 'B33'; %Euler-Bernoulli 
  %modelprops.elementtype = 'B33H' %Euler-Bernoulli
  %eltype = 'B31' %Timoshenko 
  %eltype = 'B31H' %Timoshenko 
  %modelprops.elementtype = 'B31OS'; %Timoshenko 
  %modelprops.elementtype = 'B31OSH'; %Timoshenko 
  %eltype = 'B32' %Timoshenko 
  %eltype = 'B32H' %Timoshenko 
  %modelprops.elementtype = 'B32OS'; %Timoshenko 
  %modelprops.elementtype = 'B32OSH'; %Timoshenko 
  %modelprops.elementtype = 'xx'; %current
  eltypes={'B32OSH'}
 
  
  
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
  
  modelprops.numofelm = 5;%20
  %BB5-B31OSH2048-l5-f1-eps0.001-u1
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %plotfig= [2,3,7,14,15,16,21,211,23,26,28,29]; %#ok<*NBRAK>
  %plotfig=[14,15,16,37,38,211,902];
  %plotfig=[14,15,16,900,32,7,3,43,211]; %#ok<NASGU>
  %plotfig=[15,943:945,948:954];main.savefigures=1; %#ok<NASGU>
  %plotfig=[14,15,16,952,953];
  %plotfig=[14,15,16];
  %plotfig=[2,14,35,42,47,48,50:51,53];%EV-Normierung
  plotfig=35;
  
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
  
  modelprops.epsilon = .02;  % .02;
  %epsils= {1,.5,.2,.1,.05,.02,.01,.005,.002,.001};
  %modelprops.lambda = 0:modelprops.epsilon:max(.8,30*modelprops.epsilon);%10; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.loadfactor = 1;
  %
  
  %modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=0; %-1 ... don't allow reruning, false... dont force rerun, 0.5 rerun if too less lambda, 1 force rerun
  modelprops.forcerun=0; %0 dont force, 0.5 force run if last lambda smaller than requested; 1 force run
  modelprops.numofeigs=inf;
  modelprops.allowComplex=false;
  main.closall=0;
  main.savefigures=1; % false.. dont safe figures(faster), true safe figures (slow)
  main.check=0;
  main.colorshift=5;
  modelprops.ask_delete=true;
  %modelprops.MeterValue=1; %1000mm=1m=0.001km
  main.whichEV='bungle'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; 'rNCT_K0_r';'rCT_K0_r'
  main.Normierung='A0R1'; % 'R1'; 'rCT_K0_r'; 'A0R1'
  main.rho='A0R1'; % KtR1 R1; 'A0R1'
  main.xBezug='n'; %n..normalisiert; d..differenz zut Refwert
   
  %modelprops.sigma=-10;
  modelprops.followsigma=true;
  
  % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
%   modelprops.numofelm = 256; 
%    [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%   res.stability_limit
  
  
% % % %  %% elemtype
  %eltypes={'B32OSH','B31OSH'}%,'B33','B33H','B32'}
% % %eltypes={'B32OS','B32OSH'}%,'B31OS','B31OSH'
% %eltypes={'B31OSH'}
% % %plotfig=[];
%epsils={.05,.02,.01}
epsils={modelprops.epsilon}
for j=1:numel(epsils)
    modelprops.epsilon = cell2mat(epsils(j));
    modelprops.lambda = 0:modelprops.epsilon:max(.8,30*modelprops.epsilon)
 for i=1:numel(eltypes)
  elementtype = char(eltypes(i))
  %  if strcmp(elementtype,'B32OSH') ||  strcmp(elementtype,'B31OSH')
  %   main.check=true;
  %  else
  %   main.check=false;
  %  end
  % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
  main.colorshift=(i-1)+(j-1)*numel(eltypes);
  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,[],elementtype);
  %ans=res.stability_limit
  %ans(2)*0.5e6
  res.stability_limit
 end
end
  
% %% epsilon
% epsils={.2,.1,.05,.02,.01,.005,.002}
%   for i=1:numel(epsils)
%    %modelprops.numofelm = cell2mat(i)

%    %    plotfig=[];
%    modelprops.ask_delete=true;
%    main.colorshift=i-1;
% 
% % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
% [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
% 
%   end

% % % % % % % % %numelem
% numofelms = {2,4,8,16,32,64,128,256,512};%numofelms = {2,5,10,20,50,100,200,500,1000,2000}
% for i=1:numel(numofelms)
%  modelprops.numofelm = cell2mat(numofelms(i));
%  main.colorshift=i-1;
%  
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main);
%  
%  %
% %  res.stability_limit(2)*0.5e6
%  res.stability_limit
%  
% end