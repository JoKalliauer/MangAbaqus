%#!/usr/bin/env octave
%university:TU Wien
 %#ok<*NOPTS>
 close all
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'minimized');
 set(0, 'DefaultFigureWindowState', 'normal');
 %#ok<*NBRAK> 

  % there are following predefined test cases:
  %modelprops.testcase = 'detKt2D'; % Knickstab mitte starke Achse
  %modelprops.testcase = 'detKt2Dneu'; % Knickstab mitte schwache Achse
  modelprops.testcase = 'detKt2Dgen'; % Knickstab mitte starke Achse
  %modelprops.testcase = 'd2bock';
  %modelprops.CrossSectionOrientation=[0,-1,0]%'0,0,-1\n'
  
  modelprops.length = 5;
  
  % possible element types (be aware of 2D and 3D):
  %2D:
  %eltype = 'B23'; %Euler-Bernoulli 
  %eltype = 'B23H'; %Euler-Bernoulli 
  %eltype = 'B21'; %Timoshenko
  %eltype = 'B21H'; %Timoshenko
  eltype = 'B22'; %Timoshenko 
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
  %eltypes={'B22OS','B22OSH','B21OS','B22','B21','B23','B22H'}
 
  
  
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=0;-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
 
  
  epsil =1;  % 1
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  %sortType = 'forwardJK';
  %plotfig= [15,900,913,914,908,902,916,211,7,9021,9022,906]; %#ok<*NBRAK>
  plotfig=[902]; 
  forcedeig = []; %1; % forced eigenvector number 'none' sorting
 
 
  modelprops.elementtype = eltype;
  
  %modelprops.lambda = 5*epsil; % do not go over snap-through point
  modelprops.lambda = 0:epsil:100; %(0.78-4*epsil); % do not go over snap-through point 5*epsil:10*epsil:(0.78-4*epsil)
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1;%0 if both sides
  modelprops.MeterValue=1;
  %
  
%   modelprops.a=.15;
  
  
  modelprops.forceAbaqus=0; %-1..don't run simulation even if not existing; 0 run if not existing; 1 force run
  modelprops.forcerun=0;%0..run if not exist; 0.5 run if too short; 1 force run
  modelprops.numofeigs=1;
  modelprops.allowComplex=true;
  %main.closall=true;
  main.closall=false;
  main.savefigures=uint8(1);
  main.check=true;
  %main.check=false;
  main.colorshift=0;
  modelprops.ask_delete=1;
  
  %modelprops.sigma=-5;
  modelprops.numelFac=[1 1];
  
    %modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
 %[res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,1);
%      modelprops.ask_delete=true; modelprops.forceAbaqus=false; modelprops.forcerun=false;
  
  
% list=[1,2,3,4,5,6,10,20,50];
list=[1:16,32,64,128,256,512,1024,2048];
%list=[4,8,16];
%list=[8];
%list=[5,50]
%list=[8,32,128];

% % %plotfig=[];
% 
for i=1:numel(list)
 numofelm=list(i);
 main.colorshift=i-1;
   [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,numofelm);
end


% eltypes={'B22','B21','B23','B22H','B21H','B23H'}
% % %plotfig=[];
% for i=1:numel(eltypes)
%  elementtype = char(eltypes(i))
% %  if strcmp(elementtype,'B32OSH') ||  strcmp(elementtype,'B31OSH')
% %   main.check=true;
% %  else
% %   main.check=false;
% %  end
%  % % modelprops.ask_delete=false; modelprops.forceAbaqus=true; modelprops.forcerun=true;
%  main.colorshift=i-1;
%  [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,8,[],elementtype);
%  
% end

%modelprops=rmfield(modelprops,'forceAbaqus');
