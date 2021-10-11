function [model] = runEigenProblem(modelprops)
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
 else
  NotSyncedFolder='~/Abaqus/MangAbaqus/';% NotSyncedFolder='./' 
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
 lambda0   = sort(unique(round([0,lambda]   *100000))/100000);
 
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
 fulllambda= sort(unique(round(fulllambda*100000))/100000);
 
 

 
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
   warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
   return
  end
 end
 model = selectModel(modelprops,AbaqusRunsFolder);
 model.fulllambda=fulllambda;

 %% Check if *.mat exists
 %filename = model.filename;
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem check mat-file');end
 files=dir([AnalysisResultsFolder,model.filename,'-*.mat']);
 if exist([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'], 'file') == 2 && modelprops.forcerun<=.501
  tmp=modelprops.numofeigs;
  mpl=modelprops.lambda;
  if strcmp(modelprops.testcase,'TL_arch3D')
   mpl=min(mpl,.8);
  end
  if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem load mat-file');end
  loadFileName=[AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'] %#ok<NOPRT>
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
 [Kts,num,~,BC,model.inDOF]  = AbaqusModelsGeneration.getStiffnessMatrices(model,[],modelprops.typeofanalysis);
 
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
 [~, Nres, Kg] = AbaqusModelsGeneration.getHistoryOutputFromDatFile([model.AbaqusRunsFolder,model.filename,'.dat']);
 % [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
 
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
  warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),fulllambda(min(size(Kts,1),numel(fulllambda))),max(fulllambda))
  matches(matches+2>mintest)=[];
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
 
 model.lambda0 = lambda0';
 %modelDisp=model;
 model.BC=sort(BC);
 model.fullload0=[0;model.fullload];
 model = runEigenProblemSub(modelprops,model,Displ,Kts,Kg,matches,wbrEP);
 %modelDisp=model;
 %[Kts3]  = AbaqusModelsGeneration.getStiffnessMatrices3(model,[],modelprops.typeofanalysis);
%  [~, Nres3] = AbaqusModelsGeneration.getHistoryOutputFromDatFile([model.AbaqusRunsFolder,model.filename,'.dat']);
 [Displ3,Rot3] = NodalResults2DisplJK(Nres);
 model=runEigenProblemDispJK(modelprops,model,Displ3,[],[],matches,wbrEP);
 model=runEigenProblemRotJK([],model,Rot3,[],[],matches,wbrEP);
 
 if usejava('jvm'); waitbar(1,wbrEP,'runEigenProblem EigvalueProblem finsih');end
 
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
 
 
 clear tmp ml mpl
 disp(['ready to save: ','AnalysisResults/',model.filename,'-',modelprops.typeofanalysis,'-',num2str(model.numofeigs),'.mat']);
 if isunix && ~exist(AnalysisResultsFolder, 'dir')
  mkdir(AnalysisResultsFolder)
 end
 dt=whos('model');
 if numel(model.lambda)>10
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

% function [eq1,eq2] = solcheck(Kt, typeofanal, dKt, ddKt, lam, epsil, Ls, rs, i, model)
%     eigpo = 1;
%
%     L1 = Ls(5,eigpo) + lam;
%     r01 = rs(4,:,eigpo); r01 = r01(:);  r01(model.BC(:,1)) = [];
%     r0  = rs(5,:,eigpo); r0 = r0(:);    r0(model.BC(:,1)) = [];
%     r11 = rs(6,:,eigpo); r11 = r11(:);  r11(model.BC(:,1)) = [];
%     r12 = rs(7,:,eigpo); r12 = r12(:);  r12(model.BC(:,1)) = [];
%     r13 = rs(8,:,eigpo); r13 = r13(:);  r13(model.BC(:,1)) = [];
%     r1 = r0;
%     if i == 1 % forward difference method
%         dL1 = 1/epsil*(-3/2*Ls(5,eigpo) + 2*Ls(6,eigpo) - 1/2*Ls(7,eigpo));
%         ddL1 = 1/epsil^2*(2*Ls(5,eigpo) - 5*Ls(6,eigpo) + 4*Ls(7,eigpo) - 1*Ls(8,eigpo));
%         dr1 = 1/epsil*(-3/2*r0 + 2*r11 - 1/2*r12);
%         ddr1 = 1/epsil^2*(2*r0 - 5*r11 + 4*r12 - 1*r13);
%     else % central difference method
%         dL1 = 1/epsil*(-1/2*Ls(4,eigpo) + 1/2*Ls(6,eigpo));
%         ddL1 = 1/epsil^2*(1/2*Ls(4,eigpo) -2*Ls(5,eigpo)  + 1/2*Ls(6,eigpo));
%         dr1 = 1/epsil*(-1/2*r01 + 1/2*r11);
%         ddr1 = 1/epsil^2*(1/2*r01 -2*r0  + 1/2*r11);
%     end
%
%
%     switch typeofanal
%         case 'I'
%             B = -eye(size(Kt));
%             dB = 0;
%         case 'CLE'
%             B = -dKt;
%             dB = -ddKt;
%         case 'I-CLE'
%             B = eye(size(Kt)) - dKt;
%             dB = -ddKt;
%     end
%
%     eq1   = norm((full(Kt) - (L1-lam)*B)*r1);
%     eq2_1 = norm((dKt - (dL1-1)*B - (L1-lam)*dB)*r1);
%     eq2_2 = norm((Kt - (L1-lam)*B)*dr1);
%     eq2   = norm((dKt - (dL1-1)*B - (L1-lam)*dB)*r1 + (Kt - (L1-lam)*B)*dr1);
%
% end