function model = runEigenProblem(modelprops)
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
 
<<<<<<< HEAD
 
 if usejava('jvm'); wbrEP=waitbar(0,'runEigenProblem start','Name','runEigenProblem_waitbar'); end
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 lambda =   modelprops.lambda;
 if sum(strcmp(fieldnames(modelprops), 'forcerun')) == 0
  modelprops.forcerun=true;
 end
 if sum(strcmp(fieldnames(modelprops), 'forceAbaqus')) == 0
  modelprops.forceAbaqus=false;
 elseif modelprops.forceAbaqus==true
  %modelprops.forcerun=true;
<<<<<<< HEAD
 end
 if strcmp(modelprops.testcase,'pureBendingBeamMalendowski')
  NotSyncedFolder='./';
 else
  NotSyncedFolder='~/Abaqus/MangAbaqus/';% NotSyncedFolder='./' 
 end

 if sum(strcmp(fieldnames(modelprops), 'followsigma')) == 0
  modelprops.followsigma=false;
 end
 
=======
 end
 NotSyncedFolder='~/Abaqus/MangAbaqus/';% NotSyncedFolder='./' 
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
 %lambda14 = lambda0 + 4*epsil;
 %lambda15 = []; %lambda0 + 5*epsil;
 lambda21 = max(lambda - epsil,0);
 lambda22 = max(lambda - 2*epsil,0);
 lambda23 = max(lambda - 3*epsil,0);
 %lambda24 = []; %max(lambda - 4*epsil,0);
 %lambda25 = []; %lambda - 5*epsil;
 fulllambda= [lambda0,lambda11,lambda21,lambda12,lambda22,lambda13,lambda23]';
=======
 lambda14 = lambda0 + 4*epsil;
 lambda15 = []; %lambda0 + 5*epsil;
 lambda21 = max(lambda - epsil,0);
 lambda22 = max(lambda - 2*epsil,0);
 lambda23 = max(lambda - 3*epsil,0);
 lambda24 = []; %max(lambda - 4*epsil,0);
 lambda25 = []; %lambda - 5*epsil;
 fulllambda= [lambda0,lambda11,lambda21,lambda12,lambda22,lambda13,lambda23,lambda14,lambda24,lambda15,lambda25]';
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 fulllambda= sort(unique(round(fulllambda*100000))/100000);
 
 

<<<<<<< HEAD
 
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
=======
 
 if numel(fulllambda)>405
  warning('MyProgram:Input','using %f>404 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
  %assert(numel(fulllambda)<=1005,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
%  assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
  assert(numel(fulllambda)<=37205,'using %f>2000 fulllambda-values will take more than 12GB by Abaqus',numel(fulllambda))
 end
 
 modelprops.lambda = fulllambda;
 
 
 %% creates the inputfile
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 model = selectModel(modelprops,AbaqusRunsFolder);
 model.fulllambda=fulllambda;

 %% Check if *.mat exists
 %filename = model.filename;
<<<<<<< HEAD
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem check mat-file');end
 files=dir([AnalysisResultsFolder,model.filename,'-*.mat']);
 if exist([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'], 'file') == 2 && modelprops.forcerun==false
  tmp=modelprops.numofeigs;
  %mpl=modelprops.lambda;
  if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem load mat-file');end
=======
 if exist([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'], 'file') == 2 && modelprops.forcerun==false
  tmp=modelprops.numofeigs;
  %mpl=modelprops.lambda;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  load([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'],'model');
  %only use the result if numofeigs is the same, 
  %since numofeigs changes the results
  %ml=model.lambda;
  if usejava('jvm'); close(wbrEP);end
  if tmp==size(model.eigenvalues{1},2) %&& mpl(end)==ml(end)
   clear tmp ml mpl
   return
  else
   warning('MyProgram:Input','requested %f eigenvalues, but %f exist in mat-file',tmp,size(model.eigenvalues{1},2))
  end
  clear tmp ml mpl
 elseif ~isempty(files) && modelprops.forcerun==false
  if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem search mat-file');end
  modeldef=model;
  load([AnalysisResultsFolder,files(1).name],'model');
  modelload=model;
  model=modeldef;
  nrldef=numel(model.lambda);
  nrlload=numel(modelload.lambda);
  if nrldef<nrlload
   warning('MyProgram:Inputchange','Less Lambda requested than in another run')
   model.lambda(nrldef+1:nrlload)=modelload.lambda(nrldef+1:nrlload);
   lambda0   = sort(unique(round([0,transpose(model.lambda)]   *100000))/100000);
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
<<<<<<< HEAD
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem run Abaqus');end
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 if ~exist([model.AbaqusRunsFolder,model.filename,'_STIF9.mtx'],'file')
  noresults=true;
 else
  noresults=false;
 end
 
 if (modelprops.forcerun==true && modelprops.forceAbaqus==true) || ~usejava('desktop') || noresults==true
  if usejava('desktop')
   %assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
   assert(numel(fulllambda)<=10005,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
  else
<<<<<<< HEAD
   %assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
=======
   assert(numel(fulllambda)<=2195,'using %f>504 fulllambda-values will take more than 10min by Abaqus',numel(fulllambda))
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
   assert(numel(fulllambda)<=37205,'using %f>2000 fulllambda-values will take more than 12GB by Abaqus',numel(fulllambda))
  end
  AbaqusModelsGeneration.runAbaqus(model.filename,AbaqusRunsFolder,modelprops);
 end
 % Kg = AbaqusModelsGeneration.getKgmatrix(model);
 
 %% get Stiffness
<<<<<<< HEAD
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem get Stiffness');end
 %[StifMatrices,num0,activeNodes,activeDofs,unactiveDofs,BC,Loads]  = AbaqusModelsGeneration.getStiffnessMatrices(model)
 [Kts,num,~,BC]  = AbaqusModelsGeneration.getStiffnessMatrices(model,[],modelprops.typeofanalysis);
%  if size(Kts,1)>numel(fulllambda)
%   warning('MyProgram:OldFiles','Old *.mtx-files ignored');
%   Ktsnew=Kts{1:numel(fulllambda),:};
%  end
=======
 %[StifMatrices,num0,activeNodes,activeDofs,unactiveDofs,BC,Loads]  = AbaqusModelsGeneration.getStiffnessMatrices(model)
 [Kts,num,~,BC]  = AbaqusModelsGeneration.getStiffnessMatrices(model);
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 model.BC = BC;
 
%  for j = 1:size(BC,1)
%   activeDofs(activeDofs==BC(j,1)) = []; %unactiveDofs = [unactiveDofs; BC(j,1)];
%  end
%  %unactiveDofs = sort(unique(unactiveDofs));
<<<<<<< HEAD
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem getHistoryOutputFromDatFile');end
=======
 
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 [~, Nres, EigRes] = AbaqusModelsGeneration.getHistoryOutputFromDatFile([model.AbaqusRunsFolder,model.filename,'.dat']);
 % [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
 
 Displ = NodalResults2Displ(Nres);
%  if numel(Displ)<1
%   displacementsenable=false;
%  else
%   displacementsenable=true;
%  end
 Kg = EigRes;
 
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
 
<<<<<<< HEAD
 %Energy = zeros(length(lambda0),3);
 %Energy(:,1) = lambda0;
 mintest=max(min(size(Kts,1),size(Displ,1)),3);
 if max(matches)+2>mintest
  warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),fulllambda(min(size(Kts,1),numel(fulllambda))),max(fulllambda))
  matches(matches+2>mintest)=[];
 end
=======
 Energy = zeros(length(lambda0),3);
 Energy(:,1) = lambda0;
 if max(matches)+3>size(Kts,1)
  warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),fulllambda(size(Kts,1)),max(fulllambda))
  matches(matches+3>size(Kts,1))=[];
 end
 
 %% solve EigvalueProblem
 %numofeigs=min([modelprops.numofeigs,size(Kt0,1),numofeigs0,numofeigs11,numofeigs12,numofeigs13,numofeigs14])
 %fullEV=NaN(modelprops.numofeigs,size(fulllambda,1));
 for i = 1:length(matches)
  %disp('Lambda:');
  %disp(lambda(matches(i)));
  %Lambda=fulllambda(matches(i)) %#ok<NASGU,NOPRT>
  %disp('-------------');
  %     Energy(i,2) = membrane(matches(i));
  %     Energy(i,3) = nonmembrane(matches(i));
  if i==1
   %[~,~,numofeigs] = solveCLEforMinEigNew();
   
   % 0.
   Kt0 = Kts{matches(i),2};
   %numofeigs=min(numofeigs,size(Kt0,1));
   %rr = zeros(size(Kt0,1),1);
   
   %        ru = [];
   %        ru = BC(:,1);
   %    ruMalendowski=diag(Kt0==1e36);
   %    if any(ruMalendowski)
   %     warning('MyProgram:Malendowski','Malendowski would remove some values')
   %    end
   
   ru = diag(Kt0==1e36);% remove boundary conditions
   %sizeKt0=size(Kt0,1);
   %ra = transpose(1:sizeKt0);
   %ra(ru) = [];
   
   Kt0(ru,:) = []; Kt0(:,ru) = [];
   newsizeKt0=size(Kt0,1);
   newra = transpose(1:newsizeKt0);
   %numofeigs=min(numofeigs,size(Kt0,1));
   %RR0 = zeros(sizeKt0,numofeigs);
   Kt0_0 = Kt0;
   %        Kg(ru,:) = []; Kg(:,ru) = [];
   % 1.
   dksi11 = sqrt(Displ{matches(i)}'*Displ{matches(i)});
   Kt11 = Kts{matches(i)+1,2}; Kt11(ru,:) = []; Kt11(:,ru) = [];
   % 2.
   dksi12 = sqrt((Displ{matches(i)+1}-Displ{matches(i)})'*(Displ{matches(i)+1}-Displ{matches(i)}));
   Kt12 = Kts{matches(i)+2,2}; Kt12(ru,:) = []; Kt12(:,ru) = [];
   % 3.
   dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));
   Kt13 = Kts{matches(i)+3,2};  Kt13(ru,:) = []; Kt13(:,ru) = [];
   % 4.
   dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
   Kt14 = Kts{matches(i)+4,2}; Kt14(ru,:) = []; Kt14(:,ru) = [];
   
   dksi = [0, 0, 0, 0, dksi11,dksi12,dksi13,dksi14];
   arclengths{i} = dksi;
   
   Ktprim0 = 1/epsil*(-1/2*Kt12 + 2*Kt11 - 3/2*Kt0);
   Ktprim11 = 1/epsil/2*(Kt12 - Kt0);
   Ktprim12 = 1/epsil/2*(Kt13 - Kt11);
   Ktprim13 = 1/epsil/2*(Kt14 - Kt12);
   Ktprim14 = 1/epsil*(3/2*Kt14 - 2*Kt13 + 1/2*Kt12);
   
   dKtprim0 = 1/epsil*(2*Kt0 - 5*Kt11 + 4*Kt12 - Kt13);
   
   %numofeigs=modelprops.numofeigs;
   [r0_ ,eigval0_ ,numofeigs0 ] = solveCLEforMinEigNew(Kt0,Ktprim0,Kg,Kt0_0,modelprops.typeofanalysis,matches(i),NaN,modelprops);
   [r11_,eigval11_,numofeigs11] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+1,NaN,modelprops);
   [r12_,eigval12_,numofeigs12] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+2,NaN,modelprops);
   [r13_,eigval13_,numofeigs13] = solveCLEforMinEigNew(Kt13,Ktprim13,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+3,NaN,modelprops);
   [r14_,eigval14_,numofeigs14] = solveCLEforMinEigNew(Kt14,Ktprim14,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+4,NaN,modelprops);
   
   numofeigs=min([modelprops.numofeigs,newsizeKt0,numofeigs0,numofeigs11,numofeigs12,numofeigs13,numofeigs14]);
   fullEV=NaN(numofeigs,size(fulllambda,1));
   RR0 = NaN(newsizeKt0,numofeigs);
   
   displacements_(:,5) = 0*Displ{matches(i)};
   displacements_(:,6) = Displ{matches(i)-0};
   displacements_(:,7) = Displ{matches(i)+1};
   displacements_(:,8) = Displ{matches(i)+2};
   displacements_(:,9) = Displ{matches(i)+3};
   displacements{i} = displacements_;
   
   if size(RR0,1)~=length(r0_)
    r0t = RR0; r11t = RR0; r12t = RR0; r13t = RR0; r14t = RR0;
    r0t(newra,:) = r0_;
    r11t(newra,:) = r11_;
    r12t(newra,:) = r12_;
    r13t(newra,:) = r13_;
    r14t(newra,:) = r14_;
   else
    r0t = r0_/norm(r0_);
    r11t = r11_/norm(r11_);
    r12t = r12_/norm(r12_);
    r13t = r13_/norm(r13_);
    r14t = r14_/norm(r14_);
   end
   
   
   for k = 1:size(r0_,2)
    if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
    if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
    if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end
    if r13t(:,k)'*r14t(:,k)<0; r14t(:,k) = -r14t(:,k); end
   end
   
   R = NaN(9,size(r0t,1),size(r0t,2));
   R(5,:,:) = r0t;
   R(6,:,:) = r11t;
   R(7,:,:) = r12t;
   R(8,:,:) = r13t;
   R(9,:,:) = r14t;
   if max(abs(imag(R)))>eps(0)
    warning('MyProgram:Complex','R is komplex');
   end
   
   EV = zeros(9,length(eigval0_));
   EV(5,:) = eigval0_;
   EV(6,:) = eigval11_;
   EV(7,:) = eigval12_;
   EV(8,:) = eigval13_;
   EV(9,:) = eigval14_;
   %assert(isreal(EV),'EV is complex')
   
  else % if i~= 1
   % -4.
   Zeile=matches(i)-4;
   if Zeile>0
    if Zeile>1
     dksi04 = sqrt((Displ{Zeile}-Displ{matches(i)-5})'*(Displ{Zeile}-Displ{matches(i)-5}));
     displacements_(:,1) = Displ{matches(i)-5};
    else
     dksi04 = NaN;
     displacements_(:,1) = NaN*Displ{1};
    end
    Kt04 = Kts{Zeile,2};
    dksi03 = sqrt((Displ{matches(i)-3}-Displ{matches(i)-4})'*(Displ{matches(i)-3}-Displ{matches(i)-4}));
    displacements_(:,2) = Displ{matches(i)-4};
   else
    Kt04 = 0*Kts{1,2};
    dksi04 = NaN;
    dksi03 = NaN;
    displacements_(:,1) = NaN*Displ{1};
    displacements_(:,2) = NaN*Displ{1};
   end
   Kt04(ru,:) = []; Kt04(:,ru) = [];
   % -3.
   if Zeile>-1
    Kt03 = Kts{matches(i)-3,2};
    dksi02 = sqrt((Displ{matches(i)-2}-Displ{matches(i)-3})'*(Displ{matches(i)-2}-Displ{matches(i)-3}));
    displacements_(:,3) = Displ{matches(i)-3};
   else
    Kt03 = 0*Kts{1,2};
    dksi02 = NaN;
    displacements_(:,3) = NaN;
   end
   Kt03(ru,:) = []; Kt03(:,ru) = [];
   % -2.
   if Zeile>-2
    Kt02 = Kts{matches(i)-2,2};
    dksi01 = sqrt((Displ{matches(i)-1}-Displ{matches(i)-2})'*(Displ{matches(i)-1}-Displ{matches(i)-2}));
    displacements_(:,4) = Displ{matches(i)-2};
   else
    Kt02 = 0*Kts{1,2};
    dksi01 = NaN;
    displacements_(:,4) = NaN;
   end
   Kt02(ru,:) = []; Kt02(:,ru) = [];
   % -1.
   Kt01 = Kts{matches(i)-1,2}; Kt01(ru,:) = []; Kt01(:,ru) = [];
   % 0.
   Kt0 = Kts{matches(i),2}; Kt0(ru,:) = []; Kt0(:,ru) = [];
   % 1.
   dksi11 = sqrt((Displ{matches(i)}-Displ{matches(i)-1})'*(Displ{matches(i)}-Displ{matches(i)-1}));
   Kt11 = Kts{matches(i)+1,2};
   if numel(Kt11)>0
    Kt11(ru,:) = []; Kt11(:,ru) = [];
   else
    Kt11=0*Kt0;
   end
   % 2.
   dksi12 = sqrt((Displ{matches(i)+1}-Displ{matches(i)})'*(Displ{matches(i)+1}-Displ{matches(i)}));
   Kt12 = Kts{matches(i)+2,2};
   if numel(Kt12)>0
    Kt12(ru,:) = []; Kt12(:,ru) = [];
   else
    Kt12=0*Kt11;
   end
   % 3.
   dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));
   Kt13 = Kts{matches(i)+3,2};
   if numel(Kt13)>0
    Kt13(ru,:) = [];
    Kt13(:,ru) = [];
   end
   % 4.
   
   
   
   
   
   
   
   
   displacements_(:,5) = Displ{matches(i)-1};
   displacements_(:,6) = Displ{matches(i)-0};
   displacements_(:,7) = Displ{matches(i)+1};
   displacements_(:,8) = Displ{matches(i)+2};
   
   
   
   Ktprim03 = 1/(2*epsil)*(Kt02 - Kt04);
   Ktprim02 = 1/(2*epsil)*(Kt01 - Kt03);
   Ktprim01 = 1/(2*epsil)*(Kt0 - Kt02);
   Ktprim0 = 1/(2*epsil)*(Kt11 - Kt01);
   Ktprim11 = 1/(2*epsil)*(Kt12 - Kt0);
   if numel(Kt13)>0
    Ktprim12 = 1/(2*epsil)*(Kt13 - Kt11);
   end
   
   if strcmp(modelprops.typeofanalysis,'KNL2') || size(Kts,1)<matches(i)+4
    %Kt14=0*Kt13;
    Ktprim13=0*Ktprim12;
    dksi14=NaN;
    displacements_(:,9)=NaN*displacements_(:,8);
   else
    Kt14 = Kts{matches(i)+4,2};
    Kt14(ru,:) = [];
    Kt14(:,ru) = [];
    Ktprim13 = 1/(2*epsil)*(Kt14 - Kt12);
    dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
    displacements_(:,9) = Displ{matches(i)+3};
   end
   
   dksi = [dksi04, dksi03, dksi02, dksi01, dksi11,dksi12,dksi13,dksi14];
   arclengths{i} = dksi;
   
   displacements{i} = displacements_;
   
   
   dKtprim0 = 1/epsil*(Kt01 -2*Kt0 +Kt11);
   
   [r03_,eigval03_] = solveCLEforMinEigNew(Kt03,Ktprim03,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-3,NaN,modelprops);
   [r02_,eigval02_] = solveCLEforMinEigNew(Kt02,Ktprim02,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-2,NaN,modelprops);
   [r01_,eigval01_] = solveCLEforMinEigNew(Kt01,Ktprim01,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-1,NaN,modelprops);
   [r0_ ,eigval0_ ] = solveCLEforMinEigNew(Kt0 ,Ktprim0 ,Kg,Kt0_0,modelprops.typeofanalysis,matches(i),model,modelprops); % Kt=Kt0; iter=matches(i);
   [r11_,eigval11_] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+1,NaN,modelprops);
   [r12_,eigval12_] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+2,NaN,modelprops);
   if numel(Kt13)>0
    [r13_,eigval13_] = solveCLEforMinEigNew(Kt13,Ktprim13,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+3,NaN,modelprops);
   else
    r13_=NaN*r12_;
    eigval13_=NaN*eigval12_;
   end
   
   
   
   if size(RR0,1)~=length(r0_)
    r03t = RR0; r02t = RR0; r01t = RR0; r0t = RR0; r11t = RR0; r12t = RR0; r13t = RR0;
    r03t(newra,:) = r03_;
    r02t(newra,:) = r02_;
    r01t(newra,:) = r01_;
    r0t(newra,:) = r0_;
    r11t(newra,:) = r11_;
    r12t(newra,:) = r12_;
    r13t(newra,:) = r13_;
   else
    
    r03t = r03_;%/norm(r03_);
    r02t = r02_;%/norm(r02_);
    r01t = r01_;%/norm(r01_);
    r0t = r0_;%/norm(r0_);
    r11t = r11_;%/norm(r11_);
    r12t = r12_;%/norm(r12_);
    r13t = r13_;%/norm(r13_);
   end
   
   
   for k = 1:size(r0_,2)
    if r03t(:,k)'*r02t(:,k)<0; r02t(:,k) = -r02t(:,k); end
    if r02t(:,k)'*r01t(:,k)<0; r01t(:,k) = -r01t(:,k); end
    if r01t(:,k)'*r0t(:,k)<0;  r0t(:,k)  = -r0t(:,k);  end
    if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
    if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
    if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end
   end
   
   R = NaN(9,size(r0t,1),size(r0t,2));
   R(2,:,:) = r03t;
   R(3,:,:) = r02t;
   R(4,:,:) = r01t;
   R(5,:,:) = r0t;
   R(6,:,:) = r11t;
   R(7,:,:) = r12t;
   R(8,:,:) = r13t;
   if max(abs(imag(R)))>eps(0)
    warning('MyProgram:Complex','R is komplex');
   end
   
   EV = NaN(9,length(eigval0_));
   EV(2,:) = eigval03_;
   EV(3,:) = eigval02_;
   EV(4,:) = eigval01_;
   EV(5,:) = eigval0_;
   EV(6,:) = eigval11_;
   EV(7,:) = eigval12_;
   EV(8,:) = eigval13_;
   

   
  
  end
  
  if matches(i)<4
   %first=max(matches(i)-3,1);
   fullEV(:,matches(i):matches(i)+3)=[eigval0_,eigval11_,eigval12_,eigval13_];
   %fullLAMDA(matches(i):matches(i)+3)=fulllambda(matches(i):matches(i)+3);
   if matches(i)>=2
    fullEV(:,matches(i)-1)=eigval01_;
    %fullLAMDA(matches(i):matches(i)+3)
    if matches(i)>=3
     fullEV(:,matches(i)-2)=eigval02_;
    end
   end
  else
   fullEV(:,matches(i)-3:matches(i)+3)=[eigval03_,eigval02_,eigval01_,eigval0_,eigval11_,eigval12_,eigval13_];
   %fullLAMDA(matches(i)-3:matches(i)+3)=matches(i)-3:matches(i)+3;
  end % matches(i)<4

  eigval{i} = EV;
  eigvec{i} = R;
  StiffMtxs{i,1} = Kt0;
  StiffMtxs{i,2} = Ktprim0;
  StiffMtxs{i,3} = dKtprim0;
  
  %    [eq1,eq2] = solcheck(Kt0, typeofanal, Ktprim0, dKtprim0, lambda(matches(i)), epsil, EV, R, i, model);
  %     for p = 1:10
  %         r = R(5,:,p); r = r(:);
  %         r(ru) = [];
  %         L = EV(5,p);
  %         disp(norm((Kt0 + L*Ktprim0)*r));
  %     end
 end%for i = 1:length(matches)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
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
 
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblem');end
 
 model.lambda0 = lambda0';
<<<<<<< HEAD
 model = runEigenProblemSub(modelprops,model,Displ,Kts,Kg,matches,wbrEP);
 
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
=======
 model.stiffnessMatrices = StiffMtxs;
 %model.fullLAMDA=fullLAMDA;
 model.fullEV=fullEV;
 
 clear tmp ml mpl
 disp(['ready to save: ','AnalysisResults/',model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat']);
 if isunix && ~exist(AnalysisResultsFolder, 'dir')
  mkdir(AnalysisResultsFolder)
 end
 save([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'],'model');
 dt=whos('model');
 if dt.bytes>2*1024^3
  warning('MyProgram:Size','model needs %f GB (> 2GB) space',dt.bytes/1024^3)
  save([AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'],'model','-v7.3');
 end
 disp('saved');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
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