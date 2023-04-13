%#!
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2023)
%author of subprograms: Johannes Kalliauer(2020-2023), Michal Malendowski (2019-2020)
%created: ~2020

%% Last changes
%2023-03-13 JK: main.Normierung='R1';
%2023-03-16 JK: removed old comments
%2023-04-13 JK: modelprops.allowComplex was defined twice, leading to ignoring of the first definition

%% Code

 %#ok<*NOPTS>
 %clear
 close all
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
  
 
  [~,modelprops.ecc]=eccfromU(0.5);
  modelprops.testcase = 'eccenCompressionBeam'; 

  
  modelprops.length = 5;
  
  eltypes={'B31','B31H','B32','B32H','B32OS','B32OSH'};
 

  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0;
 
  sortType = 'none';  %sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards','forwardJK'
  plotfig=[3,6,14,15,45,59];
  forcedeig = []; %1; % forced eigenvector number 'none' sorting

  
  modelprops.loadfactor = 1;
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=-1; %-1..returns error if not exist, 0..use old if exist, 1.. force new calc
  modelprops.forcerun=false; %0..use existing one, 0.5.. force run if last lambda smaller than requested, always fore a new calc.
  modelprops.allowComplex=2;%0..no complex, 1 also complex, 2 only complex
  main.closall=true;
  main.savefigures=true;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.rsame=NaN;%0.8;
  main.rstabil=NaN;%0.99999;
  main.whichEV='bungle'; % main.whichEV='bungle'; 'Disp'; 'Rot'; 'wrap'; 'Hyb'; 'bungle_rKr'; 'skip' ; 'bungle_rK0r'; 'bungle_K0r1';'rNCT_K0_r';'rCT_K0_r'; 'k11' 'k0_11'
  main.Normierung='R1'; % 'R1'; 'rCT_K0_r' 'sqrtK_r' 'skip' 'k0_11'
  main.rho='R1'; % KtR1 R1 'A0R1' 
  
  %modelprops.MeterValue=1; %1000mm=1m=0.001km ;
  main.xBezug='1'; %n..normalisiert; d..differenz zut Refwert; 1...Abaqus-Lambda; s...Stepnumber; i..individual
  main.flipAxis=false;
  
  modelprops.sigma=-1;
  modelprops.followsigma=true;
  modelprops.sortJKeigval=-1; %1..closest to zero, -1 ..most negative one
  
numofelms={2,20};


Exz={modelprops.ecc};modelprops.numofeigs=14;%min 7 EV

epsils={.02,0.005}%

for l=1:numel(epsils)
 modelprops.epsilon = cell2mat(epsils(l));
 modelprops.lambda = 0:modelprops.epsilon:max([2,20*modelprops.epsilon])
 for k=1:numel(Exz)
  for j=1:numel(numofelms)
   modelprops.numofelm = cell2mat(numofelms(j));
   for i=1:numel(eltypes)
    elementtype = char(eltypes(i))
    main.colorshift=(i-1)+(j-1)*numel(eltypes)+(l-1)*numel(eltypes)*numel(numofelms);
    [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,cell2mat(Exz(k)),elementtype);
    
   end
  end
 end
end

allFigures = findall(0,'Type','figure'); % find all figures
set(allFigures,'WindowState','normal'); % set the WindowState of all figures to normal
