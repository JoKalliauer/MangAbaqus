function [model] = runEigenProblem(modelprops)
%% run the the EigenProblem outer program
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)

%% Input
% modelprops ... parameters which were used in the Abaqus-run

%% Output
% model ... results ot the Eigenvector-Problem

%% Structure
% * runEigenProblem ... run the Eingenvalue-Problem
%   * selectModel .. calls a function to create the input-file 
%   * AbaqusModelsGeneration.runAbaqus ... run the input-file in Abaqus
%   * AbaqusModelsGeneration.getStiffnessMatrices ... get the stiffness-matrix from Abaqus-results
%   * AbaqusModelsGeneration.getHistoryOutputFromDatFile ... get the nodal-results from Abaus
%   * runEigenProblemSub ... Run the core of the eigenvalue-Problem and saving it into model
%     * solveCLEforMinEigNew ... get one specific eigenvalue and the eigenvector
%   * runEigenProblemDispJK ... Posprocessing the displacements

%% Recent Changes
%2023-02-16 JK: datetime("today") instead of date and added warning-idenfifyer, deleted solcheck by Malendowski
%2023-04-13 JK: changed name from model.lambda0 to model.lambdainput

%% Inputcheck

 %modelname,lambda,epsil,numofelm,typeofanal,additionalParameters
 if nargin<1
  %testcase = 'TL_arch';
  %numofelm = 20;
  %eltype = 'B21';
  modelprops.lambda = [.025,.075,.125,.175,.225,.2750,.325,.375,.425,.475,.525,.575,.625,.675,.725];
  modelprops.epsilon = 0.005;
  modelprops.typeofanalysis = 'K0';
  %loadFactor = 1.0;
  %len = [];
 else
  %testcase = modelprops.testcase;
  %numofelm = modelprops.numofelm;
  %eltype =   modelprops.elementtype;
  
  %typeofanalJKXY = modelprops.typeofanalysis;
  %loadFactor = modelprops.loadfactor;
  %len = [];
 end
 
 if usejava('jvm'); wbrEP=waitbar(0,'runEigenProblem start','name','runEigenProblem'); else; wbrEP=[]; end %,'WindowState','minimized'
 lambda =   modelprops.lambda;
 %loadFactor = modelprops.loadfactor;
 switch modelprops.testcase
  case 'TL_arch3D'
   RefLastMalendowski=83.3*10^3*10^2;
   if strcmp(modelprops.RefLast,'Eh')
    h = modelprops.profil.h; %20e-2;
    E=modelprops.profil.E; % 2e+11;
    RefLast=E*h;%
   else
    RefLast=modelprops.loadfactor*RefLastMalendowski;%Malendowski's Load
   end
   modelprops.loadfactor=RefLast/RefLastMalendowski;
   lambda=lambda/modelprops.loadfactor;
  case 'pureBendingBeamJK'
   if modelprops.numofelm==1
    if strcmp(modelprops.elementtype,'B31OSH')
     warning('MyPrgm:Input','B31OSH need at least 2 elements to calculate a stiffnessmatrix')
    end
   end
   MV=modelprops.MeterValue;
   RefLastMalendowski=0.5e6*MV*MV; % [Nm]
   modelprops.profil.E=2.1e+11/MV;
   if sum(strcmp(fieldnames(modelprops), 'orientate')) == 0
    modelprops.orientate=5;
    %error('not tested')
   end
   if strcmp(modelprops.RefLast,'ES')
    if modelprops.orientate==5
     Sy=653.6*10^-6*MV*MV*MV;
     %Sz=55.07*10^-6*MV*MV*MV;
     RefLast=modelprops.profil.E*Sy;
    end
   else
    RefLast=modelprops.loadfactor*0.5e6*MV*MV; %[Nm]
   end
   modelprops.loadfactor=RefLast/RefLastMalendowski;
   lambda=lambda/modelprops.loadfactor;
 end


 %% Code
 
 
 if sum(strcmp(fieldnames(modelprops), 'forcerun')) == 0
  modelprops.forcerun=true;
 end
 if sum(strcmp(fieldnames(modelprops), 'forceAbaqus')) == 0
  modelprops.forceAbaqus=false;
 elseif modelprops.forceAbaqus==true
  %modelprops.forcerun=true;
 end
 if strcmp(modelprops.testcase,'pureBendingBeamMalendowski')
  NotSyncedFolder='./';
  warning('MyProgram:outdated','pureBendingBeamMalendowski should only be used for historical reasons')
 else
   if isunix %linux-pC
    NotSyncedFolder='~/Abaqus/MangAbaqus/';% NotSyncedFolder='./' 
   elseif ispc  %windows
     warning('MyProgram:OS','This program is written for Linux, some features might not work on Windows')
     NotSyncedFolder='C:\Data\Abaqus\MangAbaqus\';% NotSyncedFolder='./'
     %NotSyncedFolder='%USERPROFILE%/Abaqus/MangAbaqus/';% NotSyncedFolder='./'
   elseif ismac %Mac-PC
     warning('MyProgram:OS','not tested on Mac-pc')
   else
     warning('MyProgram:OS','OS unknown')
   end
 end

 if sum(strcmp(fieldnames(modelprops), 'followsigma')) == 0
  modelprops.followsigma=false;
 end
 
 AbaqusRunsFolder=[NotSyncedFolder 'AbaqusRuns/'];% AbaqusRunsFolder='./AbaqusRuns/';
 %AbaqusRunsFolder='./AbaqusRuns/';
 AnalysisResultsFolder=[NotSyncedFolder 'AnalysisResults/'];
 
 %if isempty(len)
 %len = 0;
 %end
 
 lambda   = reshape(lambda,1,length(lambda));
 RundenGenauigkeit=1e8;
 lambda0   = sort(unique(round([0,lambda]   *RundenGenauigkeit))/RundenGenauigkeit);
 
 diff=min(lambda0(2:end)-lambda0(1:end-1));
 if diff+1/100000<modelprops.epsilon
  warning('MyProgram:Input','lamda-steps smaler than epsilon')
  modelprops.epsilon=min(diff,modelprops.epsilon);
 end
 epsil =    modelprops.epsilon;
 
 
 
 lambda11 = lambda0 + epsil;
 lambda12 = lambda0 + 2*epsil;
 lambda13 = lambda0 + 3*epsil;
 %lambda14 = lambda0 + 4*epsil;
 %lambda15 = []; %lambda0 + 5*epsil;
 lambda21 = max(lambda - epsil,0);
 lambda22 = max(lambda - 2*epsil,0);
 lambda23 = max(lambda - 3*epsil,0);
 %lambda24 = []; %max(lambda - 4*epsil,0);
 %lambda25 = []; %lambda - 5*epsil;
 fulllambda= [lambda0,lambda11,lambda21,lambda12,lambda22,lambda13,lambda23]';
 fulllambda= sort(unique(round(fulllambda*RundenGenauigkeit))/RundenGenauigkeit);
 
 

 
 if numel(fulllambda)>405
  warning('MyProgram:Input','using %f>404 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
  %assert(numel(fulllambda)<=1005,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
%  assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
  assert(numel(fulllambda)<=37205,'using %f>2000 fulllambda-values will take more than 12GB by Abaqus',numel(fulllambda))
 end
 
 modelprops.lambda = fulllambda;
 
 
 %% creates the inputfile
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem creates model'); end
 if ~exist(AbaqusRunsFolder, 'dir')
  if isunix
   warning('MyProgram:OS','AbaqusRunsFolder does not exist')
   mkdir(AbaqusRunsFolder);
  end
  if ispc
   warning('MyProgram:OS','AbaqusRunsFolder does not exist, try creating it on at %s',AbaqusRunsFolder)
   mkdir(AbaqusRunsFolder);
   warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder did not exist, this program was written for Linux')
   %return
  end
 end
 model = selectModel(modelprops,AbaqusRunsFolder);
 model.fulllambda=fulllambda;

 %% Check if *.mat exists
 %filename = model.filename;
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem check mat-file');end
 files=dir([AnalysisResultsFolder,model.filename,'-*.mat']);
 loadFileName=[AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'] %#ok<NOPRT>
 if     exist(loadFileName, 'file') == 2 && modelprops.forcerun<=.501
  tmp=modelprops.numofeigs;
  mpl=modelprops.lambda;
  if strcmp(modelprops.testcase,'TL_arch3D')
   mpl=min(mpl,.8);
  end
  if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem load mat-file');end
  FileInfo = dir(loadFileName);
  TimeStamp = FileInfo.date %#ok<NASGU,NOPRT>
  load(loadFileName,'model');
  %only use the result if numofeigs is the same, 
  %since numofeigs changes the results
  ml=model.lambda;
  if tmp==size(model.eigenvalues{1},2)
   if mpl(end-3)<=ml(end) || modelprops.forcerun<.499
    if usejava('jvm'); close(wbrEP);end
    clear tmp ml mpl
    return
   else
    warning('MyPrgm:Different','modelprops.lambda(end) > model.lambda(end), reruning simulation')
    if modelprops.forceAbaqus>-1
     modelprops.forceAbaqus=true;
    end
   end
  else
   warning('MyProgram:Input','requested %f eigenvalues, but %f exist in mat-file',tmp,size(model.eigenvalues{1},2))
  end
  clear tmp ml mpl
 elseif ~isempty(files) && modelprops.forcerun==false % es existiert eine ähnliche Datei aber nicht die gleiche
  if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem search mat-file');end
  modeldef=model;
  retryload=true;
  load([AnalysisResultsFolder,files(1).name],'model');
  modelload=model;
  model=modeldef;
  nrldef=numel(model.lambda);
  nrlload=min(numel(modelload.lambda),numel(modelload.load));
  if nrldef<nrlload
   warning('MyProgram:Inputchange','Less Lambda requested than in another run')
   model.lambda(nrldef+1:nrlload)=modelload.lambda(nrldef+1:nrlload);
   %lambda0   = sort(unique(round([0,transpose(model.lambda)]   *100000))/100000);
   if isfield(modelload,'load')
    model.load(nrldef+1:nrlload)=modelload.load(nrldef+1:nrlload);
   end
  end
  nrfldef=numel(model.fulllambda);
  nrflload=numel(modelload.fulllambda);
  if nrfldef<nrflload
   warning('MyProgram:Inputchange','Less fulllambda requested than in another run')
   model.fulllambda(nrfldef+1:nrflload)=modelload.fulllambda(nrfldef+1:nrflload);
   fulllambda(nrfldef+1:nrflload)=model.fulllambda(nrfldef+1:nrflload);
   
  end
 end
 
 %% Run Abaqus
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem run Abaqus');end
 if ~exist([model.AbaqusRunsFolder,model.filename,'_STIF9.mtx'],'file')
  noresults=true;
 else
  noresults=false;
 end
 

disp(['run: ','AnalysisResults/',model.filename,'-',num2str(modelprops.numofeigs)]);
  if (modelprops.forcerun>=0.499 && modelprops.forceAbaqus==true) || ~usejava('desktop') || noresults==true % wenn (a) es erzwungen wird (b) es im Terminal läuft oder (c) wenn es keine Ergebnisse gibt
   if usejava('desktop')
    %assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
    assert(numel(fulllambda)<=10005,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
   else
    %assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
    assert(numel(fulllambda)<=37205,'using %f>2000 fulllambda-values will take more than 12GB by Abaqus',numel(fulllambda))
   end
   if modelprops.forceAbaqus>=0
    AbaqusModelsGeneration.runAbaqus(model.filename,AbaqusRunsFolder,modelprops);
   else
    error('MyProgram:Input','No Abaqus results found')
   end
  end
 % Kg = AbaqusModelsGeneration.getKgmatrix(model);
 
 %% get Stiffness
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem get Stiffness');end
 %[StifMatrices,num0,activeNodes,activeDofs,unactiveDofs,BC,Loads]  = AbaqusModelsGeneration.getStiffnessMatrices(model)
 [Kts,num,~,BC,model.inDOF,model.dofs]  = AbaqusModelsGeneration.getStiffnessMatrices(model,[],modelprops.typeofanalysis,modelprops.elementtype);
 
 if exist('retryload','var')
  if retryload==true
   Kt0_0 = Kts{1,2};
   ru = diag(Kt0_0==1e36);% remove boundary conditions
   Kt0_0(ru,:) = []; Kt0_0(:,ru) = [];
   newsizeKt0=size(Kt0_0,1);
   if modelprops.numofeigs>newsizeKt0 
    model.numofeigs=min(newsizeKt0,modelprops.numofeigs);
    if exist([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat'], 'file') == 2 && modelprops.forcerun==false
     load([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat'],'model')
     warning('MyProgram:NumEigs','too many eigs requested, no further check')
     return
    end
   end
  end
 end
 
%  if size(Kts,1)>numel(fulllambda)
%   warning('MyProgram:OldFiles','Old *.mtx-files ignored');
%   Ktsnew=Kts{1:numel(fulllambda),:};
%  end
%  model.BC = BC;
 
%  for j = 1:size(BC,1)
%   activeDofs(activeDofs==BC(j,1)) = []; %unactiveDofs = [unactiveDofs; BC(j,1)];
%  end
%  %unactiveDofs = sort(unique(unactiveDofs));
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem getHistoryOutputFromDatFile');end
 [ELres, Nres, Kg] = AbaqusModelsGeneration.getHistoryOutputFromDatFile([model.AbaqusRunsFolder,model.filename,'.dat']);
 if exist('ELres','var')
  [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
  model.Energyratio=(nonmembrane)./(membrane+nonmembrane);
  model.EnergyBending=nonmembrane;
  model.EnergyMembrane=membrane;
 end
 
 Displ = NodalResults2Displ(Nres);
%  if numel(Displ)<1
%   displacementsenable=false;
%  else
%   displacementsenable=true;
%  end
 %Kg = EigRes;
 
 %eigval = cell(length(lambda0),1);
 %eigvec = cell(length(lambda0),1);
 %displacements = cell(length(lambda0),1);
 %arclengths = cell(length(lambda0),1);
 %StiffMtxs = cell(length(lambda0),3);
 %DetKtx = NaN(length(lambda0),1);
 
 matches = NaN(0);
 n = 0;
 %m = 0;
 iend=min((size(num,1)-1),size(fulllambda,1));
 for i = 1:iend
  if ~isempty(find(lambda0==fulllambda(i), 1))
   n = n+1;
   matches(n) = i;
  end
 end
 if matches(1)~=1
  matches = [1,matches] %#ok<NOPRT>
 end
 
 %Energy = zeros(length(lambda0),3);
 %Energy(:,1) = lambda0;
 mintest=max(min(size(Kts,1),size(Displ,1)),3);
 if max(matches)+2>mintest
  mintest2=max(size(Kts,1),3);
  if max(matches)+2>mintest2
   warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),fulllambda(min(size(Kts,1),numel(fulllambda))),max(fulllambda))
   matches(matches+2>mintest)=[];
  else
   warning('MyProgram:Abaqus','Displ missing (step %d), but Abaqus has (step %d, l=%f)',size(Displ,1),size(Kts,1),fulllambda(min(size(Kts,1),numel(fulllambda))))
   %matches(matches+2>mintest)=[];
  end
 end
 
 %% solve EigvalueProblem
 %Kt0_0 = Kts{matches(1),2};
 %ru = diag(Kt0_0==1e36);% remove boundary conditions
 %Kt0_0(ru,:) = []; Kt0_0(:,ru) = [];
 %newsizeKt0=size(Kt0_0,1);
 %[~ ,~ ,numofeigs0 ] = solveCLEforMinEigNew(Kt0_0,Kt0_0,Kg,Kt0_0,modelprops.typeofanalysis,matches(1),NaN,modelprops);
 %Kt11= Kts{matches(1)+1,2};
 %Kt11(ru,:) = []; Kt11(:,ru) = [];
 %Ktprim0_0 = 1/(epsil)*(Kt11 - Kt0_0);
 %[~ ,~ ,numofeigs0 ] = solveCLEforMinEigNew(Kt11,Ktprim0_0,Kg,Kt0_0,modelprops.typeofanalysis,matches(1)+1,model,modelprops,Ktprim0_0);
 %model.numofeigs=min([modelprops.numofeigs,newsizeKt0,numofeigs0]);
%  newra = transpose(1:newsizeKt0);
 %fullEV=NaN(model.numofeigs,size(fulllambda,1));
%  RR0 = NaN(newsizeKt0,numofeigs);
 
 %model.N0=size(Kt0_0,1);
 
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblem'); end
 
 model.lambda0input = lambda0';%current name
 %model.lambda0=model.lambda0input;%original name
 model.lambdainput=model.lambda0input;%better name
 %modelDisp=model;
 model.BC=sort(BC);
 model.fullload0=[0;model.fullload];
 model = runEigenProblemSub(modelprops,model,Displ,Kts,Kg,matches,wbrEP,AnalysisResultsFolder);
 %modelDisp=model;
 %[Kts3]  = AbaqusModelsGeneration.getStiffnessMatrices3(model,[],modelprops.typeofanalysis);
 %  [~, Nres3] = AbaqusModelsGeneration.getHistoryOutputFromDatFile([model.AbaqusRunsFolder,model.filename,'.dat']);
 if isempty(Nres.keys)
  Displ3=[];
  Rot3=[];
  warning('MyPrgm:MSG','No Displ-Results, check *.msg-file')
  %return
 else
  [Displ3,Rot3] = NodalResults2DisplJK(Nres);
 end
 if ~isempty(Displ3) && numel(Displ3)>2
  model=runEigenProblemDispJK(modelprops,model,Displ3,[],[],matches,wbrEP);
  model=runEigenProblemRotJK([],model,Rot3,[],[],matches,wbrEP);
 end
 
 if usejava('jvm'); waitbar(1,wbrEP,'runEigenProblem EigvalueProblem finish');end
 
 %model.eigenvalues = eigval;
 %model.eigenvectors = eigvec;
 %model.energy = Energy;
 %model.arclengths = arclengths;
 %model.displacements = displacements;
 %model.lambda0 = lambda0';
 %model.stiffnessMatrices = StiffMtxs;
 %model.DetKtx=DetKtx;
 %model.N=size(Kt0,2);
 assert(model.N0==model.N,'size of Kt0 changes or is not symmetric')
 %model.fullEV=fullEV;
 %model.numofeigs=numofeigs;
 model.fulllambda1=fulllambda(2:end);
 
 clear tmp ml mpl
 disp(['ready to save: ','AnalysisResults/',model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat']);
 if isunix && ~exist(AnalysisResultsFolder, 'dir')
  mkdir(AnalysisResultsFolder)
 end
 dt=whos('model');
 if numel(model.lambda)>10
  model.date=datetime("today");%model.date=date;
  save([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat'],'model');
  if dt.bytes>2*1024^3
   warning('MyProgram:Size','model needs %f GB (> 2GB) space',dt.bytes/1024^3)
   save([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat'],'model','-v7.3');
  end
 end
 disp('saved');
 if usejava('jvm'); waitbar(i/length(matches),wbrEP,'runEigenProblem saved');close(wbrEP);end
 
%  if any(all(isnan(real(model.fullEV)))) %wenn es Spalten gibt die in jeder Zeile immer NaN sind
%   model.fulllambda(all(isnan(real(model.fullEV)), 1)) = [];
%   model.fullEV(:,all(isnan(real(model.fullEV)), 1)) = [];
%  end


end %fucntion
