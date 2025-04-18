function model = runEigenProblemSub(modelprops,model,Displ,Kts,Kg,matches,wbrEP,AnalysisResultsFolder)
%#!/usr/bin/env octave -q
%university:TU Wien
%author:Johannes Kalliauer(©2020-2023), Michał Malendowski (©2019-2020)

%% run the core of the EigenProblem

%% Input
% modelprops ... parameters which were used in the Abaqus-run
% model ... structure for the Matlab-results
% Displ ... Diplacements of the beam
% Kts ... Stiffness matrix
% Kg ... only relevant for modelprops.typeofanalysis == 'Kg' (outdated)
% matches ... which lambda refers to which results
% wbrEP .. Waitbar for Eigenproblem
% AnalysisResultsFolder ... Folder of the Abaqus-Results

%% Output
% model ... results ot the Eigenvector-Problem


%% Recent Changes
%2023-02-16 JK: added comments for explanation, added error-identifyer
%2023-02-21 JK: added k0_11
%2023-04-13 JK: changed name from model.lambda0 to model.lambdainput
%2023-09-09 JK: added r_dotKB1_r

%% Code
 
% [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);

%% Input-check
%Displ = NodalResults2Displ(Nres);
if numel(Displ)<1
 displacementsenable=false;
 if isempty(Displ)
  %emptyDispl=true;
 else
  warning('Myprgm:Unexpected','unexpected behaviour')
 end
elseif numel(Displ)==1
 warning('myPrgm:Sub:VARmissing','no disp data?')
 displacementsenable=false; 
else
 displacementsenable=true;
end


%% inizializing
dofs=model.dofs([1 2]);
if strcmp(modelprops.whichEV,'k0_11') || strcmp(modelprops.whichEV,'k11') ||  strcmp(modelprops.whichEV,'2023-12')
 HybridNodes=model.inDOF(2)-model.inDOF(1)+1;
end
%Kg = EigRes;
lambda0=model.lambdainput';
lenLam0=length(lambda0);
lenMatch=length(matches);
eigval = cell(lenLam0,1);
eigvec = cell(lenLam0,1);
eigvec2023 = cell(lenLam0,1);
eigvec1 = cell(lenLam0,1);
displacements = cell(lenLam0,1);
darclengths = cell(lenLam0,1);
arclengthJK = cell(lenLam0,1);
arclengthHM = cell(lenLam0,1);
StiffMtxs = cell(lenLam0,3);
DetKtx = NaN(lenLam0,1);
%load0 = NaN(lenLam0,1);
eigvecDRH = cell(lenLam0,1);% DRH...Displacement,Rotation,Hybrid(splitted)
eigvecDR= cell(lenLam0,1);
eigvecH= cell(lenLam0,1);
eigvecH2 = cell(lenLam0,1);% DRH...Displacement,Rotation,Hybrid(splitted)



rddotKtr = NaN(lenMatch,1);
ZweirKtt = NaN(lenMatch,1);
%NormR1=NaN(lenMatch,1);
rdotKtr = NaN(lenMatch,1);
rdotKtt = NaN(lenMatch,1);
rKtr = NaN(lenMatch,1);
rKt0r = NaN(lenMatch,1);
RerKt0Imr = NaN(lenMatch,1);
NormKt0r = NaN(lenMatch,1);
KB1Klammer = cell(lenMatch,1);
dotKB1 = cell(lenMatch,1);
ddotKB1 = cell(lenMatch,1);
t_KB1_t = NaN(lenMatch,1);
imagValues= NaN(lenMatch,1);
RrA0Rr = NaN(lenMatch,1);
RrARr = NaN(lenMatch,1);
r_dotKB1_r = NaN(lenMatch,1);
r_dotKB1_t = NaN(lenMatch,1);
r_ddotKB1_r = NaN(lenMatch,1);

%  matches = NaN(0);
%  n = 0;
%m = 0;
%  iend=min((size(num,1)-1),size(fulllambda,1));
%  for i = 1:iend
%   if ~isempty(find(lambda0==fulllambda(i), 1))
%    n = n+1;
%    matches(n) = i;
%   end
%  end
%  if matches(1)~=1
%   matches = [1,matches] %#ok<NOPRT>
%  end

%Energy = zeros(lenLam0,3);
%Energy(:,1) = lambda0;
mintest=max(min(size(Kts,1),size(Displ,1)),3);
%  if max(matches)+2>mintest
%   %warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),model.fulllambda(min(size(Kts,1),numel(model.fulllambda))),max(model.fulllambda))
%   matches(matches+2>mintest)=[];
%  end

%% solve EigvalueProblem
Kt0_0 = Kts{matches(1),2};
oldsizeKt0=size(Kts{matches(1),2},1);
ru36 = diag(Kt0_0==1e36);% remove boundary conditions
ruinf = diag(Kt0_0==inf);
ru=logical(ru36+ruinf);
BC=find(ru);
aktiveDOF=dofs(1)*dofs(2)-numel(BC);% which dofs are displacements and Rotations
assert(all(BC==model.BC(:,1)),'ru does not agree with BC')
Kt0_0(ru,:) = []; Kt0_0(:,ru) = [];
newsizeKt0=size(Kt0_0,1);
Kt11= Kts{matches(1)+1,2};
Kt11(ru,:) = []; Kt11(:,ru) = [];
Ktprim0_0 = 1/(modelprops.epsilon)*(Kt11 - Kt0_0);
[r0t0 ,~ ,numofeigs0 ] = solveCLEforMinEigNew(Kt11,Ktprim0_0,Kg,Kt0_0,modelprops.typeofanalysis,matches(1)+1,model,modelprops,Ktprim0_0);
numofeigs=min([modelprops.numofeigs,newsizeKt0]);
if numofeigs<numofeigs0
 warning('MyPrgm:unexpexted','numofeigs<numofeigs0')
end
if numofeigs<modelprops.numofeigs
 warning('MyPrgm:NumEigs','less eigs available')
 if modelprops.forcerun==0
  modelprops.numofeigs=numofeigs;
  loadFileName=[AnalysisResultsFolder,model.filename,'-',modelprops.typeofanalysis,'-',num2str(modelprops.numofeigs),'.mat'] %#ok<NOPRT>
  if     exist(loadFileName, 'file') == 2
   load(loadFileName,'model');
   return
  end
 end
end
newra = transpose(1:newsizeKt0);
fullEV=NaN(numofeigs,size(model.fulllambda,1));
RR0 = NaN(newsizeKt0,numofeigs);
DHtmp=NaN(newsizeKt0,numofeigs);
Dtmp=NaN(aktiveDOF,numofeigs);
%t=NaN(newsizeKt0,numofeigs);
tPunkt01=NaN(newsizeKt0,numofeigs);
tPunkt=NaN(newsizeKt0,numofeigs);
tPunkt11=NaN(newsizeKt0,numofeigs);
%eigvecA0r = cell(lenLam0,numofeigs);
NormeigvecA0r = NaN(lenLam0,numofeigs);
NormR1= NaN(lenLam0,numofeigs);
rKt0rij= NaN(lenLam0,numofeigs);
rCTKt0rij= NaN(lenLam0,numofeigs);%https://de.mathworks.com/help/matlab/ref/transpose.html
rNCTKt0rij= NaN(lenLam0,numofeigs);%https://de.mathworks.com/help/matlab/ref/ctranspose.html
RerNCTKt0Rerij= NaN(lenLam0,numofeigs);
RerCTKt0Rerij= NaN(lenLam0,numofeigs);
RerARerij= NaN(lenLam0,numofeigs);
cosphirij = NaN(lenLam0,numofeigs);
if numel(Displ)>0
 sizeDisp=size(Displ{1},1);
 if size(Displ{1},1)~=oldsizeKt0
  if ~strcmp(modelprops.elementtype(end),'H')%not shure what this warning means
   warning('MyProgram:Disp','size(Displ{1},1) has a different size than oldsizeKt0')
  end
 end
else
 warning('MyProgram:Disp','Displ from *.dat are empty, maybe Abaqus-Error')
 if ~strcmp(modelprops.testcase,'d2bockDisp')
  assert(0,'Displ are empty, check *.msg-file')
 end
 sizeDisp=min(oldsizeKt0);
end
displacements_=NaN(sizeDisp,9);

model.N0=size(Kt0_0,1);

if strcmp(modelprops.elementtype(1:2),'B2')
 model.RestrictedDOFs=3:5;
 model.Dim=2;
else
 model.RestrictedDOFs=[];
 model.Dim=3;
end

R = NaN(7,size(r0t0,1),size(r0t0,2)); % dl x DoF x NrEigs
R_DRsize=[5,model.dofpNode,size(model.Nodes,1),size(r0t0,2)]; % dl x DoFpNode x Nodes x NrEigs
HNodes=model.inDOF(2)-model.inDOF(1)+1;
%R_Hsize=[5,3,HNodes,size(r0t0,2)];
R_DRHsize=R_DRsize+[0 -numel(model.RestrictedDOFs) HNodes 0]; % dl x DoFpNode(reducedRestricted) x Nodes(inkl.Hyb) x NrEigs
R_DRH = NaN(R_DRHsize); % increments x DoFpNode x Nodes x NrEigs
R_DRH2023 = NaN(R_DRHsize([2 3 1 4])); % DoFpNode x Nodes x increments x NrEigs
R_DRHsizeIn=R_DRHsize(2:4);%+[0 inNr 0]; %DoFpNode x Nodes(inkl.Hyb) x NrEigs
ReducedHybridDofa=model.inDOF(3)+1:R_DRHsizeIn(1);
ReducedHybridDofb=model.inDOF(4)+1:R_DRHsizeIn(1);
nA=false(R_DRHsizeIn(1:2));
nA(ReducedHybridDofa,model.inDOF(1):2:model.inDOF(2))=true;
if ~isnan(ReducedHybridDofb)
 nA(ReducedHybridDofb,model.inDOF(1)+1:2:model.inDOF(2))=true;
end
if sum(strcmp(fieldnames(model), 'JC')) == 0
 model.JC=[];
else
 for i1=1:size(model.JC,1)
  nA(model.JC(i1,2),model.JC(i1,1))=true;%only works if model.JC(i,2) is smaller than model.RestrictedDOFs
 end
end
if sum(strcmp(fieldnames(modelprops), 'sigma')) == 0
 modelprops.sigma=0;
end

if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblem');end

taktuell=5;
f=length(matches);
if strcmp(modelprops.elementtype,'B33') || strcmp(modelprops.elementtype,'B33H') || strcmp(modelprops.elementtype,'B31H')
 isB32=false;
else
 isB32=true;
end

dotKB1{1}=NaN*Kt0_0;
ddotKB1{1}=NaN*Kt0_0;

for ilambda = 1:f 
 if usejava('jvm'); waitbar(ilambda/length(matches),wbrEP,'runEigenProblem EigvalueProblem');end
 %disp('Lambda:');
 %disp(lambda(matches(i)));
 %Lambda=fulllambda(matches(i)) %#ok<NASGU,NOPRT>
 %disp('-------------');
 %     Energy(i,2) = membrane(matches(i));
 %     Energy(i,3) = nonmembrane(matches(i));
 %for i
 %if i>0
 % -4.
 Zeile=matches(ilambda)-4;
 if Zeile>0 && displacementsenable
  if Zeile>1
   dksi04 = sqrt((Displ{Zeile}-Displ{matches(ilambda)-5})'*(Displ{Zeile}-Displ{matches(ilambda)-5}));
   displacements_(:,1) = Displ{matches(ilambda)-5};
  else
   dksi04 = NaN;
   displacements_(:,1) = NaN*Displ{1};
  end
  %Kt04 = Kts{Zeile,2};
  dksi03 = sqrt((Displ{matches(ilambda)-3}-Displ{matches(ilambda)-4})'*(Displ{matches(ilambda)-3}-Displ{matches(ilambda)-4}));
  displacements_(:,2) = Displ{matches(ilambda)-4};
 else
  %Kt04 = 0*Kts{1,2};
  dksi04 = NaN;
  dksi03 = NaN;
  if displacementsenable
   displacements_(:,1) = NaN*Displ{1};
   displacements_(:,2) = NaN*Displ{1};
  end
 end
 %Kt04(ru,:) = []; Kt04(:,ru) = [];
 % -3.
 if Zeile>-1 && displacementsenable
  Kt03 = Kts{matches(ilambda)-3,2};
  dksi02 = sqrt((Displ{matches(ilambda)-2}-Displ{matches(ilambda)-3})'*(Displ{matches(ilambda)-2}-Displ{matches(ilambda)-3}));
  displacements_(:,3) = Displ{matches(ilambda)-3};
 else
  Kt03 = 0*Kts{1,2};
  dksi02 = NaN;
  displacements_(:,3) = NaN;
 end
 Kt03(ru,:) = []; Kt03(:,ru) = [];
 % -2.
 if Zeile>-2 && displacementsenable
  Kt02 = Kts{matches(ilambda)-2,2};
  dksi01 = sqrt((Displ{matches(ilambda)-1}-Displ{matches(ilambda)-2})'*(Displ{matches(ilambda)-1}-Displ{matches(ilambda)-2}));
  displacements_(:,4) = Displ{matches(ilambda)-2};
 else
  Kt02 = 0*Kts{1,2};
  dksi01 = NaN;
  displacements_(:,4) = NaN;
 end
 Kt02(ru,:) = []; Kt02(:,ru) = [];
 % -1.
 if Zeile>-3 && displacementsenable
  Kt01 = Kts{matches(ilambda)-1,2};
  dksi11 = sqrt((Displ{matches(ilambda)}-Displ{matches(ilambda)-1})'*(Displ{matches(ilambda)}-Displ{matches(ilambda)-1}));
  displacements_(:,5) = Displ{matches(ilambda)-1};
 else
  Kt01 = 0*Kts{1,2};
  dksi11 = NaN;
  displacements_(:,5) = NaN;
 end
 Kt01(ru,:) = []; Kt01(:,ru) = [];
 % 0.
 KT = Kts{matches(ilambda),2}; KT(ru,:) = []; KT(:,ru) = [];
 % 1.
 Kt11 = Kts{matches(ilambda)+1,2};
 if numel(Kt11)>0
  Kt11(ru,:) = []; Kt11(:,ru) = [];
 else
  Kt11=0*KT;
 end
 % 2.
 if displacementsenable
  dksi12 = sqrt((Displ{matches(ilambda)+1}-Displ{matches(ilambda)})'*(Displ{matches(ilambda)+1}-Displ{matches(ilambda)}));
 else
  dksi12 = NaN;
 end

 % 3.
 %dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));

 % 4.
 if displacementsenable
  displacements_(:,6) = Displ{matches(ilambda)-0};
  displacements_(:,7) = Displ{matches(ilambda)+1};
  if size(Displ,1)<matches(ilambda)+2
   displacements_(:,8)=NaN*displacements_(:,7);
  else
   displacements_(:,8) = Displ{matches(ilambda)+2};
  end
  %displacements_(:,9) = Displ{matches(i)+3};
 end

 %Ktprim03 = 1/(2*epsil)*(Kt02 - Kt04);
 Ktprim02 = 1/(2*modelprops.epsilon)*(Kt01 - Kt03);
 Ktprim01 = 1/(2*modelprops.epsilon)*(KT - Kt02);
 if any(Kt01(:))
  Ktprim0 = 1/(2*modelprops.epsilon)*(Kt11 - Kt01);
 else
  Ktprim0 = 1/(modelprops.epsilon)*(Kt11 - KT);
 end
 Kt2prim0 = 1/(modelprops.epsilon^2)*(Kt11+Kt01-2*KT);



 if size(Kts,1)<matches(ilambda)+2
  Ktprim11=0*Ktprim0;
 else
  Kt12 = Kts{matches(ilambda)+2,2};
  if numel(Kt12)>0
   Kt12(ru,:) = []; Kt12(:,ru) = [];
  else
   Kt12=0*Kt11;
  end
  Ktprim11 = 1/(2*modelprops.epsilon)*(Kt12 - KT);
 end

 if size(Kts,1)<matches(ilambda)+3
  Ktprim12=0*Ktprim11;
 else
  Kt13 = Kts{matches(ilambda)+3,2};
  if numel(Kt13)>0
   Kt13(ru,:) = [];
   Kt13(:,ru) = [];
   Ktprim12 = 1/(2*modelprops.epsilon)*(Kt13 - Kt11);
  else
   Ktprim12=0*Ktprim11;
  end
 end

 if strcmp(modelprops.typeofanalysis,'KNL2') || mintest<matches(ilambda)+4
  %Kt14=0*Kt13;
  %Ktprim13=0*Ktprim12;
  %dksi14=NaN;
  if displacementsenable; displacements_(:,9)=NaN*displacements_(:,8);end
 else
  %Kt14 = Kts{matches(i)+4,2};
  %Kt14(ru,:) = [];
  %Kt14(:,ru) = [];
  %Ktprim13 = 1/(2*epsil)*(Kt14 - Kt12);
  %Ktprim14 = 1/epsil*(3/2*Kt14 - 2*Kt13 + 1/2*Kt12);
  %dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
  displacements_(:,9) = Displ{matches(ilambda)+3};
 end

 dksi = [dksi04, dksi03, dksi02, dksi01, dksi11,dksi12];
 darclengths{ilambda} = dksi;

 displacements_(abs(displacements_)<=1e-31)=0;%remove numeric issues close to zero
 displacements{ilambda} = displacements_;
 arclengthJK{ilambda}=mean(sqrt(displacements_(:,5).*displacements_(:,5)));
 arclengthHM{ilambda}=sqrt(sum(displacements_(:,5).*displacements_(:,5)));

 if ilambda>1
  if matches(ilambda)==matches(ilambda-1)+1
   foloworder=true;
  else
   foloworder=false;
  end
 else
  foloworder=false;
 end

 %dKtprim0 = 1/epsil*(Kt01 -2*Kt0 +Kt11);
 if strcmp(modelprops.typeofanalysis,'Kg')
  error('MyPrgm:removed','this function got removed, use older version')
 end

  if foloworder==true
   %r03=r02_
   r02_=r01_;
   r01_=r0_;
   r0_=r11_;
   r11_=r12_;
   eigval02_=eigval01_;
   eigval01_=eigval0_;
   eigval0_=eigval11_;
   eigval11_=eigval12_;
  else
   %[r03_,eigval03_] = solveCLEforMinEigNew(Kt03,Ktprim03,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-3,NaN,modelprops);
   [r02_,eigval02_] = solveCLEforMinEigNew(Kt02,Ktprim02,Kg,Kt0_0,modelprops.typeofanalysis,matches(ilambda)-2,model,modelprops,Ktprim0_0);
   [r01_,eigval01_] = solveCLEforMinEigNew(Kt01,Ktprim01,Kg,Kt0_0,modelprops.typeofanalysis,matches(ilambda)-1,model,modelprops,Ktprim0_0);
   [r0_ ,eigval0_,~,KB1Klammer{ilambda},imagValues(ilambda)] = solveCLEforMinEigNew(KT ,Ktprim0 ,Kg,Kt0_0,modelprops.typeofanalysis,matches(ilambda),model,modelprops,Ktprim0_0); % Kt=Kt0_0; iter=matches(i); typeofanal=modelprops.typeofanalysis;
   [r11_,eigval11_,~,KB1Klammer{ilambda+1},imagValues(ilambda+1)] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,modelprops.typeofanalysis,matches(ilambda)+1,model,modelprops,Ktprim0_0);
  end
  if ~any(Ktprim12(:)) && strcmp(modelprops.typeofanalysis,'CLE')
   error('MyPrgm:removed','this function got removed, use older version')
  end
   if ~any(Kt12) % Check if there are still values to process
    if ilambda~=length(matches)
     warning('MyPrgm:Zero','Kt12 is zero')
    end
    r12_=NaN*r11_;
    eigval12_=NaN*eigval11_;
   else %there are still values to process
    [r12_,eigval12_,~,KB1Klammer{ilambda+2},imagValues(ilambda+2)] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,modelprops.typeofanalysis,matches(ilambda)+2,model,modelprops,Ktprim0_0);
    %if any(Kt01(:)) && i>1
     dotKB1{ilambda+1}=1/(2*modelprops.epsilon)*(KB1Klammer{ilambda+2}-KB1Klammer{ilambda});
     ddotKB1{ilambda+1}=1/(modelprops.epsilon^2)*(KB1Klammer{ilambda+2}-2*KB1Klammer{ilambda+1}+KB1Klammer{ilambda});
    %end
   end
  
  



  if size(RR0,1)~=length(r0_)
   %r03t = RR0;
   r02t = RR0; r01t = RR0; r0t = RR0; r11t = RR0; r12t = RR0; %r13t = RR0; % r14t = RR0;
   %r03t(newra,:) = r03_;
   r02t(newra,:) = r02_;
   r01t(newra,:) = r01_;
   r0t(newra,:) = r0_;
   r11t(newra,:) = r11_;
   r12t(newra,:) = r12_;
   %r13t(newra,:) = r13_;
   %r14t(newra,:) = r14_;
  else
   %r03t = r03_;%/norm(r03_);
   r02t = r02_;%/norm(r02_);
   r01t = r01_;%/norm(r01_);
   r0t = r0_;%/norm(r0_);
   r11t = r11_;%/norm(r11_);
   r12t = r12_;%/norm(r12_);
   %r13t = r13_;%/norm(r13_);
   %r14t = r14_/norm(r14_);
  end


  for k = 1:size(r0_,2)
   if r01t(:,k)'*r0t(:,k)<0;  r01t(:,k)  = -r01t(:,k);  end
   if r02t(:,k)'*r01t(:,k)<0; r02t(:,k) = -r02t(:,k); end
   %if r03t(:,k)'*r02t(:,k)<0; r03t(:,k) = -r03t(:,k); end
   if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
   if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
   %if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end
   %if r13t(:,k)'*r14t(:,k)<0; r14t(:,k) = -r14t(:,k); end
  end

  %R(2,:,:) = r03t;
  R(taktuell-2,:,:) = r02t;
  R(taktuell-1,:,:) = r01t;
  R(taktuell,:,:) = r0t;
  R(taktuell+1,:,:) = r11t;
  R(taktuell+2,:,:) = r12t;
  %R(8,:,:) = r13t;
  %R(9,:,:) = r14t;
  if isB32
   R_DRH(1,:,:,:)=RealEV(r02t,BC,R_DRHsize,nA);
   R_DRH(2,:,:,:)=RealEV(r01t,BC,R_DRHsize,nA);
   R_DRH(3,:,:,:)=RealEV(r0t ,BC,R_DRHsize,nA);
   R_DRH(4,:,:,:)=RealEV(r11t,BC,R_DRHsize,nA);
   R_DRH(5,:,:,:)=RealEV(r12t,BC,R_DRHsize,nA);
   R_DRH2023(:,:,1,:)=RealEV(r02t,BC,R_DRHsize,nA);
   R_DRH2023(:,:,2,:)=RealEV(r01t,BC,R_DRHsize,nA);
   R_DRH2023(:,:,3,:)=RealEV(r0t ,BC,R_DRHsize,nA);
   R_DRH2023(:,:,4,:)=RealEV(r11t,BC,R_DRHsize,nA);
   R_DRH2023(:,:,5,:)=RealEV(r12t,BC,R_DRHsize,nA);
  else
   %keeping them NaN (not implemented)
   %R_DRH=NaN(5,?,?,?);
   %R_DRH2023=NaN(DoFpNode,Nodes,5,NrEigs); %DoFpNode x Nodes x increments x NrEigs
  end
  if max(abs(imag(R)))>eps(0)
   warning('MyProgram:Complex','R is komplex');
  end

  EV = NaN(7,length(eigval0_)); %overwrite everything with NaNs
  %EV(2,:) = eigval03_;
  EV(3,:) = eigval02_;
  EV(4,:) = eigval01_;
  EV(taktuell,:) = eigval0_;
  EV(6,:) = eigval11_;
  EV(7,:) = eigval12_;
  %EV(8,:) = eigval13_;
  %EV(9,:) = eigval14_;
  %EV
 if ilambda==12
  disp(r0t)
 end

 if strcmp(modelprops.typeofanalysis,'KNL2') && modelprops.followsigma==true
  modelprops.sigma=min(modelprops.sigma,min(real(EV(:))));
 end
 if modelprops.followsigma==true && numofeigs==1
  if modelprops.numofeigs==1
   modelprops.sigma=min(modelprops.sigma,min(real(eigval0_)));
  elseif ~any(isnan(eigval0_))
   modelprops.sigma=min(real(eigval0_));
  end
 end



 %end %if i>0

 if exist('eigval11_','var')
  if matches(ilambda)<4
   %first=max(matches(i)-3,1);
   fullEV(:,matches(ilambda):matches(ilambda)+2)=[eigval0_,eigval11_,eigval12_];
   %fullLAMDA(matches(i):matches(i)+3)=fulllambda(matches(i):matches(i)+3);
   if matches(ilambda)>=2
    fullEV(:,matches(ilambda)-1)=eigval01_;
    %fullLAMDA(matches(i):matches(i)+3)
    if matches(ilambda)>=3
     fullEV(:,matches(ilambda)-2)=eigval02_;
    end
   end
  else
   fullEV(:,matches(ilambda)-2:matches(ilambda)+2)=[eigval02_,eigval01_,eigval0_,eigval11_,eigval12_];
   %fullLAMDA(matches(i)-3:matches(i)+3)=matches(i)-3:matches(i)+3;
  end % matches(i)<4
 else
  fullEV(:,matches(ilambda))=eigval0_(1);
 end


 eigval{ilambda} = EV;
 %imagValues(i)=imagValuesi;
 eigvecDRH{ilambda}=R_DRH;% DRH...[Displacement,Rotation,Hybrid](splitted)  % increments x DoFpNode x Nodes x NrEigs
 %eigvecDRH2023{i}=R_DRH2023;% DRH...[Displacement,Rotation,Hybrid](splitted)  % increments x DoFpNode x Nodes x NrEigs
 eigvecDRi=R_DRH2023(:,1:R_DRsize(3),:,:);% DRH...[Displacement,Rotation,Hybrid](splitted)
 eigvecDR{ilambda}=eigvecDRi;
 if ~strcmp(modelprops.elementtype,'B21H')
  eigvecHi=R_DRH2023(1:6,R_DRsize(3)+1:end,:,:);% DoFpNode x Nodes x increments x NrEigs
  eigvecH{ilambda}=eigvecHi;
  eigvecH2i=squeeze(sum(eigvecHi.^2,2:3)); % increments x NrEigs
  eigvecH2{ilambda}=eigvecH2i;
 end
 if strcmp(modelprops.whichEV,'k0_11') || strcmp(modelprops.whichEV,'2023-12') || strcmp(modelprops.whichEV,'2023_12half')
  for incriment=3:7% 3... lambda-2*dlambda; 4 ... previous loadstep; 5 ... current loadstep; 6 next loadstep; 7 second next loadstep
   for EVNri=1:numofeigs
    DHtmp(:,EVNri)=sqrt(full(diag(Kt0_0))).'.*R(incriment,:,EVNri);% sqrt(k_ii)*r_i
    Dtmpi=DHtmp(1:aktiveDOF,EVNri);%only Displacement&Rotation no Hybrid
    if all(isnan(Dtmpi(:)))
     Dtmpi=NaN*Dtmpi;%useless line
     if ilambda>2
      warning('MyPrgm:NaNi','ilambda=%f, Dtmp is nan',ilambda)
     end
    elseif strcmp(modelprops.whichEV,'2023_12half')
     Dtmpi=Dtmpi/norm(Dtmpi);
    end
    Dtmp(:,EVNri)=Dtmpi;
    %t(:,EVNri)=R(incriment,:,EVNri);%
    if all(isnan(R(4,:,EVNri)))
     tPunkt(:,EVNri)=NaN*(R(6,:,EVNri)-R(4,:,EVNri))/(2*modelprops.epsilon);%to prevent strange behaviour with real nan but nonNan imaginary
    else
     tPunkt01(:,EVNri)=(R(5,:,EVNri)-R(3,:,EVNri))/(2*modelprops.epsilon);
     tPunkt(:,EVNri)=(R(6,:,EVNri)-R(4,:,EVNri))/(2*modelprops.epsilon);
     tPunkt11(:,EVNri)=(R(7,:,EVNri)-R(5,:,EVNri))/(2*modelprops.epsilon);
     if strcmp(modelprops.whichEV,'2023_12half')
      tPunkt01(:,EVNri)=tPunkt01(:,EVNri)/norm(tPunkt01(:,EVNri));
      tPunkt(:,EVNri)  =  tPunkt(:,EVNri)/norm(  tPunkt(:,EVNri));
      tPunkt11(:,EVNri)=tPunkt11(:,EVNri)/norm(tPunkt11(:,EVNri));
     end %if strcmp(modelprops.whichEV,'2023_12half')
    end %if all(isnan(R(4,:,EVNri)))
   end %for EVNri=1:numofeigs
   eigvec2023{ilambda}(1:aktiveDOF,incriment-2,:)=modelprops.alphaDRW*Dtmp;%save only the displacements and Rotations (no hybrid)
  end
  if strcmp(modelprops.elementtype,'B33') || strcmp(modelprops.elementtype,'B33H')
   error('MyPrgm:Missing','elementtype not implemented')
  else
   if  strcmp(modelprops.whichEV,'k0_11')
    if all(isnan(modelprops.alphaH))
     modelprops.alphaH=norm(DHtmp(1:aktiveDOF,1))/norm(reshape(eigvecHi(:,:,4,1),[HybridNodes*6 1 1]));
    end
    eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,:,:)=modelprops.alphaH*reshape(eigvecHi(:,:,:,:),[HybridNodes*6 5 numofeigs]);% add hybrid DOFs
   elseif strcmp(modelprops.whichEV,'2023-12') || strcmp(modelprops.whichEV,'2023_12half')
    eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,2,:)=modelprops.alphaH*tPunkt01(aktiveDOF+1:newsizeKt0,:);
    eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,3,:)=modelprops.alphaH*tPunkt(aktiveDOF+1:newsizeKt0,:);
    eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,4,:)=modelprops.alphaH*tPunkt11(aktiveDOF+1:newsizeKt0,:);
    eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,1,:)=NaN*eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,2,:);
    %eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,5,:)=NaN*eigvec2023{ilambda}(aktiveDOF+1:newsizeKt0,4,:);
    if ilambda==3
     disp('debugging point')
    end 
   else
    error('MyPrgm:Missing','whichEV not implemented')
   end
  end
  if numofeigs>1
   eigvec{ilambda} = NaN*(R); % dl x DoF x NrEigs % vector not used for this normalization any more therfore saved with NaN
  else
   eigvec{ilambda} = NaN*sqrt(diag(Kt0_0)).*R(5,:)+NaN;% vector not used for this normalization any more therfore saved with NaN
  end
 elseif  strcmp(modelprops.whichEV,'g_durch_k')%möglicher name, kannst auch was anderes nehmen
  %2023-08-07: @Antonia hier würde ich die neue Normierung mit gll und sqrt(kll) vom 2023-08-07 einfügen.
 else
  eigvec{ilambda} = R; % dl x DoF x NrEigs
  if ~strcmp(modelprops.elementtype,'B21H')
   eigvec{ilambda}(:,1:aktiveDOF,:)=modelprops.alphaDRW*eigvec{ilambda}(:,1:aktiveDOF,:);
   eigvec{ilambda}(:,aktiveDOF+1:newsizeKt0,:)=modelprops.alphaH*eigvec{ilambda}(:,aktiveDOF+1:newsizeKt0,:);
  else
   warning('MyPrgm:Missing','not implemnted')
   eigvec{ilambda}=eigvec{ilambda}*NaN;
  end
 end
 StiffMtxs{ilambda,1} = KT;
 StiffMtxs{ilambda,2} = Ktprim0;
 if modelprops.numofelm<=20 % dont calculate DetKt for large, because it is slow
  if ilambda==1
   %    KTmult0=5*10^-11;
   KTmult0=1.1e-11;
   KTrel=KT*KTmult0;
   res=det(KTrel);
   relative=false;
   if (res==0 && ~strcmp(modelprops.elementtype(end),'H')) || isinf(res)
    Kt0=KT;
    Kt0inv=inv(Kt0);
    relative=true;
   else
    KTmult=KTmult0;
   end
  end %i=1
  if relative==true
   if numel(KT)<=1E5
    DetKtx(ilambda)=det(Kt0inv*KT); %#ok<MINV>
   else
    DetKtx(ilambda)=NaN;
   end
  else
   KTrel=KT*KTmult;
   if isinf(res) || res==0
    DetKtx(ilambda)=res;
   else
    DetKtx(ilambda)=det(KTrel);
   end
  end
 else
  DetKtx(ilambda)=NaN;
 end % if modelprops.numofelm<=20

 %% Normierung
 if strcmp(modelprops.whichEV,'NoHyb')
  if strcmp(modelprops.elementtype(end),'H')
   HybrideFreiheitsgrade=(diag(Kt0_0)<=0);
   r01t(HybrideFreiheitsgrade,:)=0;
   r0t(HybrideFreiheitsgrade,:)=0;
   r11t(HybrideFreiheitsgrade,:)=0;
  else
   modelprops.whichEV='bungle'; % for nonhybrid elements nothing to change
   %HybrideFreiheitsgrade=[];
  end
 end % if strcmp(modelprops.whichEV,'NoHyb')
 
 s_old={'bungle_rK0r','bungle','bungle_K0r1','NoHyb','k0_11'} ;
 s_2023={'k0_11','2023-12','2023_12half'};
 if any(strcmp(modelprops.whichEV,s_old)) 
  % if ~isfield(modelprops,'forcedeig')
  %  modelprops.forcedeig=1; % if it does not exit, create it
  % end
  if isempty(modelprops.forcedeig)
   forcedeig=1;
  else
   forcedeig=modelprops.forcedeig(1); %take a shorter abbrevation
  end
  rm=r0t(:,forcedeig);
  if strcmp(modelprops.whichEV,'bungle_rK0r') || strcmp(modelprops.whichEV,'bungle_K0r1') ||  strcmp(modelprops.whichEV,'bungle') || strcmp(modelprops.whichEV,'NoHyb') || strcmp(modelprops.whichEV,'k0_11')
   r01=r01t(:,forcedeig);
   r11=r11t(:,forcedeig);

   if strcmp(modelprops.whichEV,'bungle_rKr') || strcmp(modelprops.whichEV,'bungle_rCT_K_r') || strcmp(modelprops.whichEV,'k0_11')
    if strcmp(modelprops.whichEV,'bungle_rKr')
     warning('MyProgram:UnclearInput','please use bungle_r.K0r or bungle_rTK0r instead of bungle_rK0r')
    end
    Nenner01=sqrt(r01'*Kt01*r01);
    Nenner0=sqrt(rm'*KT*rm);
    Nenner11=sqrt(r11'*Kt11*r11);
   elseif strcmp(modelprops.whichEV,'bungle_r.K0r') %Skalarprodukt aus <r,r1>
    Nenner01=sqrt(transpose(r01)*Kt0_0*r01);
    Nenner0=sqrt(transpose(rm)*Kt0_0*rm);
    Nenner11=sqrt(transpose(r11)*Kt0_0*r11);
   elseif strcmp(modelprops.whichEV,'bungle_rTK0r') || strcmp(modelprops.whichEV,'bungle_rK0r') || strcmp(modelprops.Normierung,'rCT_K0_r') %Transponiert: r^T*r1
    if strcmp(modelprops.whichEV,'bungle_rK0r')
     warning('MyProgram:UnclearInput','please use bungle_r.K0r or bungle_rTK0r instead of bungle_rK0r')
    end
    Nenner01=sqrt(r01'*Kt0_0*r01);
    Nenner0=sqrt(rm'*Kt0_0*rm);
    Nenner11=sqrt(r11'*Kt0_0*r11);
   elseif strcmp(modelprops.whichEV,'bungle_K0r1')
    Nenner01=NaN;
    Nenner0=norm(Kt0_0*rm);
    Nenner11=NaN;
   elseif strcmp(modelprops.Normierung,'k11') || strcmp(modelprops.whichEV,'k11')
    Nenner01=sqrt(dot(r01,diag(Kt01).*r01)+eigvecH2i(2,forcedeig));
    Nenner0=sqrt(dot(rm,diag(KT).*rm)+eigvecH2i(3,forcedeig));
    Nenner11=sqrt(dot(r11',diag(Kt11).*r11)+eigvecH2i(4,forcedeig));
   elseif strcmp(modelprops.Normierung,'k0_11')
    Nenner01=sqrt(dot(r01 ,diag(Kt0_0).*r01)+eigvecH2i(2,forcedeig));
    Nenner0 =sqrt(dot(rm  ,diag(Kt0_0).*rm )+eigvecH2i(3,forcedeig));
    Nenner11=sqrt(dot(r11',diag(Kt0_0).*r11)+eigvecH2i(4,forcedeig));
   elseif strcmp(modelprops.Normierung,'R1')
    Nenner01=NaN;
    Nenner0=norm(rm);
    Nenner11=NaN;
   elseif strcmp(modelprops.Normierung,'skip')
    warning('MyPrgm:Unexpected','Due to perfomance-reasons this code should not be called')
    Nenner01=NaN;
    Nenner0=NaN;
    Nenner11=NaN;
   else
    assert(0,'not implemented')
   end

   assert(numel(Nenner11)==1,'dimension wrong');
   r01 = r01/Nenner01;
   rm = rm/Nenner0;
   r11 = r11/Nenner11;

   rKt2ri=transpose(rm)*Kt2prim0*rm;
   if isempty(rKt2ri)
    rddotKtr(ilambda)=NaN;
   else
    rddotKtr(ilambda)=rKt2ri;
   end
   v = (r11 - r01)/(2*modelprops.epsilon);
   ZweirKtt(ilambda)=transpose(rm)*KT*v+transpose(v)*KT*rm;
   %NormR1(i)=norm(rm);
   if ~isfield(modelprops,'Normierung'); warning('Myprgm:Strange','modelprops.Normierung unset assuming R1'); modelprops.Normierung='R1'; end
   if ~strcmp(modelprops.Normierung,'R1')
    assert(norm(rm)~=1,'Norm should not be one')
   end
   rdotKtr(ilambda)=transpose(rm)*Ktprim0*rm;
   rdotKtt(ilambda)=(transpose(rm)*Ktprim0*v+transpose(v)*Ktprim0*rm)/2;
   t_KB1_t(ilambda)=transpose(v)*KB1Klammer{ilambda}*v;
   RrA0Rr(ilambda)=transpose(real(rm))*Kt0_0*real(rm);%2022-05-09
   RrARr(ilambda)=transpose(real(rm))*KT*real(rm);%2022-05-09
   r_dotKB1_r(ilambda)=transpose(rm)*dotKB1{ilambda}*rm;%2023-09-09
   r_dotKB1_t(ilambda)=(transpose(rm)*dotKB1{ilambda}*v+transpose(v)*dotKB1{ilambda}*rm)/2;%2023-09-09
   r_ddotKB1_r(ilambda)=transpose(rm)*ddotKB1{ilambda}*rm;%2023-09-09
  elseif strcmp(modelprops.whichEV,'bungle')
   if ~isnan(rm)
    assert(abs(norm(rm)-1)<1e-7,'rm not 1')
   end
  else
   error('Myprgm:missing','not implemented')
  end
  rKtr(ilambda)=transpose(rm)*KT*rm;
  rKt0r(ilambda)=single(transpose(rm)*Kt0_0*rm);
  RerKt0Imr(ilambda)=single(real(transpose(rm))*Kt0_0*imag(rm));
  NormKt0r(ilambda)=single(norm(Kt0_0*rm));
 elseif strcmp(modelprops.whichEV,'sqrtK_r') || strcmp(modelprops.whichEV,'sqrtK0_r')
  %everythings fine
 elseif  any(strcmp(modelprops.whichEV,s_2023))
  r0t=NaN*r0t;
 else
  warning('MyPrgm:NotImplemented','modelprops.whichEV=%s unknown?',modelprops.whichEV)
 end%if strcmp(modelprops.whichEV,'bungle_rK0r') || strcmp(modelprops.whichEV,'bungle') || strcmp(modelprops.whichEV,'bungle_K0r1') || strcmp(modelprops.whichEV,'NoHyb')
 
 
 for j=1:numofeigs
  locNorm=strcmp(modelprops.Normierung,'rNCT_K0_r') ||  strcmp(modelprops.Normierung,'rCT_K0_r');
  if locNorm || strcmp(modelprops.whichEV,'bungle')|| strcmp(modelprops.whichEV,'split') || strcmp(modelprops.whichEV,'corrected') || strcmp(modelprops.whichEV,'sqrtK_r') || strcmp(modelprops.whichEV,'NoHyb')
   rmj=r0t(:,j);
   if ~any(isnan(rmj))
    if strcmp(modelprops.whichEV,'bungle_rKr')
     Nenner0=sqrt(rmj'*KT*rmj);
    elseif strcmp(modelprops.whichEV,'bungle_r.K0r') || strcmp(modelprops.whichEV,'rNCT_K0_r') || strcmp(modelprops.Normierung,'rNCT_K0_r') %nonconjungate Skalarprodukt aus <r,r1> might be komplex
     Nenner0=sqrt(transpose(rmj)*Kt0_0*rmj);
    elseif strcmp(modelprops.whichEV,'bungle_rTK0r') || strcmp(modelprops.whichEV,'rCT_K0_r')  ||  strcmp(modelprops.Normierung,'rCT_K0_r') %Konunjiert Komplex transponiert
     Nenner0=sqrt(ctranspose(rmj)*Kt0_0*rmj);
    elseif strcmp(modelprops.whichEV,'bungle_K0r1')  ||  strcmp(modelprops.Normierung,'A0R1')
     Nenner0=norm(Kt0_0*rmj);
    elseif strcmp(modelprops.whichEV,'sqrtK_r')
     %this seems to be useless, modify sortEigenValuesAndGetQuantities instead
     %rmj_old=rmj;
     %rmj=sqrt(full(diag(KT))).*rmj_old;
     rmj=NaN*rmj;
     Nenner0=norm(rmj);
    elseif strcmp(modelprops.Normierung,'R1') &&  ( strcmp(modelprops.whichEV,'bungle') ||  strcmp(modelprops.whichEV,'NoHyb') )
     %warning('MyPrgm:unused','this code should be unused, please check if modelprops.whichEV=%s is correct',modelprops.whichEV)
     Nenner0=norm(rmj);
    elseif strcmp(modelprops.Normierung,'k11')
     Nenner0=sqrt(dot(rm,diag(KT).*rm));
    elseif strcmp(modelprops.Normierung,'k0_11')
     Nenner0=sqrt(dot(rm,diag(Kt0_0).*rm));
    elseif strcmp(modelprops.Normierung,'skip')
     warning('MyPrgm:runEigen:Inconsistent:Input','Normierung is set to skip, but whichEV is not skip')
     warning('off','MyPrgm:runEigen:Inconsistent:Input')
     Nenner0=NaN;
    else
     assert(0,'not implemented')
    end

    if real(Nenner0)/abs(imag(Nenner0))>3e10
      Nenner0=real(Nenner0);
    end

    rmj = rmj/(Nenner0);

    eigvecA0rij=Kt0_0*(rmj);
    NormeigvecA0rij=norm(eigvecA0rij);
    NormeigvecA0r(ilambda,j)=NormeigvecA0rij;
    NormR1ij=norm(rmj);
    NormR1(ilambda,j)=NormR1ij;

    rNCTKt0rijIJ=transpose(rmj)*Kt0_0*rmj;%nonconjungate Skalar
    if real(rNCTKt0rijIJ)/abs(imag(rNCTKt0rijIJ))>7e9
     rNCTKt0rij(ilambda,j)=real(rNCTKt0rijIJ);%nonconjungate sklarprodukt
    elseif real(rNCTKt0rijIJ)<0
     warning('MyPrgm:negativ:rSKt0rijIJ','rSKt0rijIJ ist negativ')
     warning('off','MyPrgm:negativ:rSKt0rijIJ')
     rNCTKt0rij(ilambda,j)=rNCTKt0rijIJ;%nonconjungate sklarprodukt
    else
     warning('MyPrgm:Complex:rSKt0rijIJ','rSKt0rijIJ ist komplex')
     warning('off','MyPrgm:Complex:rSKt0rijIJ')
     rNCTKt0rij(ilambda,j)=NaN;%rSKt0rijIJ;%sklarprodukt
    end
    RerNCTKt0RerijIJ=transpose(real(rmj))*Kt0_0*real(rmj);%nonconjungate Sklalar
    RerNCTKt0Rerij(ilambda,j)=RerNCTKt0RerijIJ;


    rCTKt0rijIJ=ctranspose(rmj)*Kt0_0*rmj;%conjungate TRansponiert
    if real(rCTKt0rijIJ)/abs(imag(rCTKt0rijIJ))>3e10
     rCTKt0rij(ilambda,j)=real(rCTKt0rijIJ);%conjungate TRansponiert
    elseif real(rCTKt0rijIJ)<0
     warning('MyPrgm:negativ:rTKt0rijIJ','rTKt0rijIJ ist negativ')
     warning('off','MyPrgm:negativ:rTKt0rijIJ')
     rCTKt0rij(ilambda,j)=real(rCTKt0rijIJ);%conjungate TRansponiert
    else
     warning('MyPrgm:Complex:rTKt0rijIJ','rTKt0rijIJ ist komplex')
     warning('off','MyPrgm:Complex:rTKt0rijIJ')
     rCTKt0rij(ilambda,j)=NaN;%rTKt0rijIJ;%conjungate TRansponiert
    end
    RerCTKt0RerijIJ=ctranspose(real(rmj))*Kt0_0*real(rmj);%conjungate TRansponiert
    RerCTKt0Rerij(ilambda,j)=RerCTKt0RerijIJ;
    
    
    RerARerij(ilambda,j)=transpose(real(rmj))*KT*real(rmj);



    if any(transpose(rmj)~=ctranspose(rmj)) && ~any(isnan(rmj))
     warning('MyPrgm:Complex:rmj','rmj ist komplex')
     warning('off','MyPrgm:Complex:rmj')
    end

    cosphirIJij=double(rCTKt0rijIJ)/(NormR1ij*NormeigvecA0rij);
    if real(cosphirIJij)<0
     warning('MyPrgm:negativ:RealcosphirIJij','real(cosphirIJij) ist negativ')
     warning('off','MyPrgm:negativ:RealcosphirIJij')
    end
    if real(cosphirIJij)>1
     warning('MyPrgm:negativ:RealcosphirIJij','real(cosphirIJij) ist larger than 1')
    end

    cosphirij(ilambda,j)=cosphirIJij;


   else %any(isnan(rmj))
    rCTKt0rij(ilambda,j)=NaN;
    rNCTKt0rij(ilambda,j)=NaN;
    cosphirij(ilambda,j)=NaN;
   end
  elseif strcmp(modelprops.whichEV,'sqrtK0_r')
   %everything is fine, no warning, it is just implemented in sortEigenValuesAndGetQuantities and not here.
  elseif ~strcmp(modelprops.whichEV,'skip') && numofeigs>0
   %warning('MyPrgm:NotImplemented','if modelprops.whichEV=%s is implented in sortEigenValuesAndGetQuantities then please modify the code here',modelprops.whichEV)
  end
 end %j=1:numofeigs
 
 if ~strcmp(modelprops.whichEV,'skip') && numofeigs>0
  if ~exist('rmj','var')
   if ~strcmp(modelprops.whichEV,'sqrtK0_r') && ~any(strcmp(modelprops.whichEV,s_2023))
    warning('MyPrgm:Unknown','modelprops.whichEV=%s might be unkown',modelprops.whichEV)
   end
   rmj=NaN*r0t(:,j);
  end
  eigvec1{ilambda} = rmj;
 end
end%for i = 1:length(matches)
for ilambda=f+1:min(f+3,lenLam0)
StiffMtxs{ilambda}=StiffMtxs{ilambda-1}*NaN;
end

 EVsize=size(eigvec{1});
 if numel(EVsize)>2
  assert(EVsize(1)==7,'dl should be 7')
  assert(EVsize(3)==numofeigs,'numofeigs does not match')
 end
 %assert(EVsize(1)==EVsize(2),'maybe should be symmetiric?')

if usejava('jvm'); waitbar(1,wbrEP,'runEigenProblem EigvalueProblem finsih');end

% model=M1;
model.eigenvalues =( eigval);
model.imagValues=single(imagValues);
model.arclengths =( darclengths);
model.arclengthurJK =( arclengthJK);
model.arclengthurHM =( arclengthHM);
%model.displacements =( displacements);%much diskspace
model.lambda0 =single( lambda0');
model.load0=single([0;model.load]);
model.DetKtx=single(DetKtx);
model.N=single(size(KT,2));
assert(model.N0==model.N,'size of Kt0 changes or is not symmetric')
model.fullEV=single(fullEV);
model.numofeigs=single(numofeigs);
model.rddotKtr=single(rddotKtr);
model.ZweirKtt=single(ZweirKtt);
model.NormR1 =(NormR1);%single leads to jumps in function
model.rdotKtr=single(rdotKtr);
model.rdotKtt=single(rdotKtt);
model.rKtr   =single(rKtr);
model.rKt0r  =single(rKt0r);
model.RerKt0Imr=single(RerKt0Imr);
%model.KB1=KB1Klammer;%much diskspace
model.t_KB1_t=single(t_KB1_t);
model.RrA0Rr=single(RrA0Rr);
model.RrARr=single(RrARr);
model.r_dotKB1_r=single(r_dotKB1_r);
model.r_dotKB1_t=single(r_dotKB1_t);
model.r_ddotKB1_r=single(r_ddotKB1_r);
if ~strcmp(modelprops.whichEV,'skip')
 model.NormKt0r=single(NormKt0r);
 model.eigenvectors = (eigvec); % dl x DoF x NrEigs %much diskspace
 model.eigvec1=eigvec1;
 model.NormeigvecA0r=(NormeigvecA0r);%single leads to jumps in function
 model.rKt0rij=single(rKt0rij);%outdated thats NaN
 model.rCTKt0rij=single(rCTKt0rij);% model.rTKt0rij=single(rCTKt0rij);
 model.rNCTKt0rij=single(rNCTKt0rij);% model.rSKt0rij=single(rNCTKt0rij);
 model.RerNCTKt0Rerij=single(RerNCTKt0Rerij);
 model.RerCTKt0Rerij=single(RerCTKt0Rerij);
 model.RerARerij=single(RerARerij);
 model.cosphirij=(cosphirij);%single leads to jumps in function
 if sum(strcmp(fieldnames(modelprops), 'rho')) == 0
  modelprops.rho=[];
 end
 model.eigvecDR=eigvecDR;
 splitDRH = (strcmp(modelprops.whichEV,'Disp_rK0r') || strcmp(modelprops.whichEV,'Disp') || strcmp(modelprops.whichEV,'corrected') || strcmp(modelprops.whichEV,'split') || strcmp(modelprops.whichEV,'Rot')  || strcmp(modelprops.whichEV,'Hyb'));
 if strcmp(modelprops.rho,'KtR1') || strcmp(modelprops.Normierung,'rCT_K0_r') || strcmp(modelprops.Normierung,'KtR1') || strcmp(modelprops.Normierung,'A0R1') ...
   ||  strcmp(modelprops.whichEV,'sqrtK0_r')
  model.stiffnessMatrices = (StiffMtxs(1:2,1:2));%much diskspace
 elseif strcmp(modelprops.whichEV,'sqrtK_r')
  model.stiffnessMatrices = StiffMtxs;% very much diskspace
 elseif splitDRH
  model.eigvecDRH=(eigvecDRH);% DRH...Displacement,Rotation,Hybrid(splitted) %much diskspace % increments x DoFpNode x Nodes x NrEigs
 end 
 s2={'k11','k0_11','2023-12','2023_12half'};
 if any(strcmp(modelprops.whichEV,s2)) || strcmp(modelprops.Normierung,'k11') || strcmp(modelprops.Normierung,'k0_11')
  if strcmp(modelprops.whichEV,'k11')
   model.stiffnessMatrices = StiffMtxs(:,1);% very much diskspace
  else
   model.stiffnessMatrices = StiffMtxs(1:2,1);% very much diskspace
  end
  %model.eigvecDRH=(eigvecDRH);% DRH...Displacement,Rotation,Hybrid(splitted) %much diskspace % increments x DoFpNode x Nodes x NrEigs
  %model.eigvecDR=eigvecDR;%only Displacements and Rotations
  if ~strcmp(modelprops.elementtype,'B21H')
   model.eigvecH=eigvecH;%only Hybrid-DOF
   model.eigvecH2=eigvecH2;%sum of Hybrid-DOFs^2
  end
  model.eigvec2023=eigvec2023;
 else
  %model.stiffnessMatrices = StiffMtxs(1,1);
 end
end




end %fucntion

