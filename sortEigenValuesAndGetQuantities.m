function res = sortEigenValuesAndGetQuantities(model,sortType,~,forcedeig,limit,~,main) % forcedeig=k3;
%% calculate rho based on the solution of the EW-problem
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)

%% Input
% model ... results from runEigenProblem
% sortType ... how to sort the eivenvalues (not working => ignore)
% noplot,~ .. removed variable got replaced by plotfig, kept for backwards-compatibility
% forcedeig ... force a specific Eigenvalue
% limit ... if values are outside those limits consider it as numerical error
% ~ ... removed variable
% main ... parameters how to process the data

%% Output
% res ... results for rho and others


%% Recent Changes
%2023-02-16 JK: added explantation
%2023-02-21 JK: added k0_11

%% Code


if nargin<1
 model = runEigenProblem();
 sortType = 'backwards';
 %         sortType = 'forwards';
 %         sortType = 'none';
 %plotfig = [1,2];
 forcedeig = [];
end


 if strcmp(main.whichEV,'bungle_rKr') || strcmp(main.Normierung,'k11') || strcmp(main.whichEV,'k11') || strcmp(main.whichEV,'k0_11')
  StiffMtxs=model.stiffnessMatrices;% very much diskspace
 end

epsilon = model.lambda(2) - model.lambda(1);
if ~exist('limit','var')
 limit.new=true;
end
if ~exist('main','var')
%  if sum(strcmp(fieldnames(main), 'whichEV')) == 0
  main.whichEV='bungle'; % main.whichEV='bungle'; main.whichEV='Disp'; main.whichEV='Rot'; main.whichEV='wrap'; main.whichEV='Hyb'; main.whichEV='bungle_rKr';
%  end
end
if sum(strcmp(fieldnames(model), 'numofeigs')) == 0
 sizeNumEigs=size(model.eigenvalues);
 model.numofeigs=sizeNumEigs(2);
end
if sum(strcmp(fieldnames(model), 'fulllambda')) == 0
 model.fulllambda=model.lambda;
end
if sum(strcmp(fieldnames(main), 'rho')) == 0
 main.rho='R1';
end


if ~isfield(limit,'OC1')
 if epsilon<0.1
  %limit.OC1=0.0004;
  limit.OC1=0.0008;%ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
 else
  limit.OC1=0.004;
 end
end
if ~isfield(limit,'OC5')
 if epsilon<0.1
  limit.OC5=0.0002; %ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2
 %limit.OC5=4.85e-5;% ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1 epsil = 0.01 % 3.6e-5%0.00006;
 limit.OC5=abs( -0.2628); %ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
 else
  limit.OC5=0.010;% c-B32OSH-2-l5-LF-1-Steg0.0086-eps0.02
 end
end
if ~isfield(limit,'OC5a')
 if epsilon<0.05
  %limit.OC5a=0.0000044; %=0.0000044;% pureBendingBeam-B32OSH-2-len-5-loadfac--1 epsil = 0.01  % %0.000003
  %0.010*epsilon; %ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.02
  %0.2033e-03 %ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
  limit.OC5a=max(0.010*epsilon,0.2033e-03);
 elseif epsilon>=.1
  %limit.OC5a=0.0480*epsilon;
  limit.OC5a=0.0480*epsilon; %ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.1
 else
  limit.OC5a=0.0005;
 end
 limit.OC5=max(limit.OC5,limit.OC5a);
end
if ~isfield(limit,'OC6')
 if epsilon<0.1
  %old 4.2e-5
  %0.0007 from ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2
  %0.0064*epsilon from ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.02 
  limit.OC6=max(abs(-0.0007291696954872),0.0064*epsilon);
 else
  limit.OC6=abs(-0.003);
 end
end
if ~isfield(limit,'OC6a')
 limit.OC6a=5e-11;
end
if ~isfield(limit,'OC7')
 if epsilon<0.1
  %8.3e-9 old
  %-0.1482e-6 %ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
  %tmp=5.1600e-06;
  tmp=1.459e-05; %AnalysisResults/ecc-B32OS-5-l5-e0.5-f1-eps0.05-KNL2-1.mat
  limit.OC7=max(abs(-0.1482e-6),tmp*epsilon); %-5.1600e-06*epsilon from 
 else
  limit.OC7=abs(-0.59e-4);% ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-KNL2 epsil =0.1
 end
end
if ~isfield(limit,'C2')
 limit.C2=124;
end
if ~isfield(limit,'C3A0')
 %limit.C3A0=3.63;
 %limit.C3A0=5.0972;%ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
 limit.C3A0=5.4006284680542; %AnalysisResults/ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-I-2.mat
end
if ~isfield(limit,'C4At')
 %limit.C4At=3.63;
 limit.C4At=5.0604; %ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
end
if ~isfield(limit,'C4Atdiff')
 if epsilon>0.1
  limit.C4Atdiff=abs(-0.0576);
 else
  %limit.C4Atdiff=abs(-0.19);
  %limit.C4Atdiff=abs(-0.6731);%ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
  limit.C4Atdiff=+69.229134971915*epsilon;
 end
 
end
if ~isfield(limit,'C5')
 %limit.C5=2.43;
 %limit.C5=3.0176; %ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1.
 limit.C5=3.745; %AnalysisResults/ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-I-2.mat
end
if ~isfield(limit,'C5diff')

 if epsilon>=.01
  %limit.C5diff=-0.1025; % limit.C5diff=min(-0.0007,-0.07*epsilon);%-0.000017
  limit.C5diff=-0.109; %ecc-B32OS-2-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2-KNL2-1
 else
  %Feps=-1.1;% -1.1*epsilon from ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.02
  Feps=-3.332; % limit.C5diff=-0.016656 pureBendingBeam-B32OSH-20-len-5-loadfac-1-eps0.005-KNL2-1.mat 
  limit.C5diff=min(-0.0007,Feps*epsilon);  %old -0.000017
 end
else
 if limit.C5diff>0
  warning('MyProgramm:Input','limit.C5diff should be negativ')
  limit.C5diff=-inf;
 end
end
if ~isfield(limit,'C6Tau')
 %limit.C6Tau=3.9041e7;
 %limit.C6Tau=191;
 %limit.C6Tau=3000; %AnalysisResults/TL_arch3D-B32-20-loadfac-1-eps0.005-KNL2-1.mat
 %limit.C6Tau=3844;
 %limit.C6Tau=6408;%AnalysisResults/TL_arch3D-B33H-10-loadfac-1-eps0.01-KNL2-1.mat %(?)
 %limit.C6Tau=43452;%AnalysisResults/TL_arch3D-B33H-10-loadfac-1-eps0.01-KNL2-1.mat (?)
 limit.C6Tau=16339;%AnalysisResults/TL_arch3D-B32H-20-loadfac-1-eps0.005-KNL2-1.mat (?,max)
 limit.C6Tau=5757;%AnalysisResults/TL_arch3D-B32H-20-loadfac-1-eps0.01-KNL2-1.mat (?,max)
end


if ~isfield(limit,'C8minrho')
 limit.C8minrho=0;
end
limit.C9=19.6;


if sum(strcmp(fieldnames(model), 'check')) == 0
 model.check=true;
end
%model.filename
eigvalTMP=model.eigenvalues;
s2={'bungle','sqrtK_r','sqrtK0_r','NoHyb','k0_11','2023-12','2023_12Hyb','2023_12half'};
if any(strcmp(main.whichEV,s2)) ||  strcmp(main.Normierung,'rNCT_K0_r') ||  strcmp(main.Normierung,'rCT_K0_r')
 eigvecTMP = model.eigenvectors; %defined in runEigenProblemSub
 realistic=false;
 evmiddle=5;
elseif strcmp(main.whichEV,'Disp') ||  strcmp(main.whichEV,'Disp_rK0r')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 if model.Dim==2
  relDofs=1:2;
 else
  relDofs=1:3;
 end
 relNodes=1:model.inDOF(1)-1;
elseif strcmp(main.whichEV,'Rot')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 if model.Dim==2
  relDofs=3;
 else
  relDofs=4:6;
 end
 relNodes=1:model.inDOF(1)-1;
elseif strcmp(main.whichEV,'wrap')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 if model.Dim==2
  relDofs=[];
 else
  relDofs=7:model.dofpNode;
 end
 relNodes=1:model.inDOF(1)-1;
elseif strcmp(main.whichEV,'Hyb')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 relDofs=1:model.inDOF(3);
 relNodes=model.inDOF(1):model.inDOF(2);
elseif strcmp(main.whichEV,'all') || strcmp(main.whichEV,'corrected')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 relDofs=1:model.dofpNode;
 relNodes=1:model.inDOF(2);
else
%  eigvecTMP = model.eigvecDRH;
%  realistic=true;
%  evmiddle=3;
 %relDofs=...;
 %relNodes=...;
 error('MyProgam:whichEV','main.whichEV=%s not regogniced',main.whichEV)
end
EVtmpsize=size(eigvecTMP{1});
eigvec=cell(size(eigvecTMP));
DimsEVtmp=numel(EVtmpsize);
if model.numofeigs<=1.5 %Wenn man nur einen Eigenwert hat, dann schmeißt Matlab die letzte Dimension weg, daher eine Dimension mehr
 DimsEVtmp=DimsEVtmp+1;
end
if DimsEVtmp>3 %Wenn man den aufgespalteten Eigenvektor in Knoten *Freiheitsgrade nimmt
%  eigvec=cell(size(eigvecTMP));
 for ievTMP=1:numel(eigvecTMP)
  if ~isempty(eigvecTMP{ievTMP}) && exist('relDofs','var')
   if exist('relNodes','var')
    if isempty(relNodes)
     warning('MyProgram:Input','no relevant data found')
     res=[];
     return
     %forcedeig=[];norm(r0)
     %error('MyProgram:Input','no relevant data found')
    end
   else
    warning('MyPrgm:Input','relNodes undefined')
   end
   eigvec{ievTMP}=eigvecTMP{ievTMP}(:,relDofs,relNodes,:);
  else
   eigvec{ievTMP}=NaN*eigvecTMP{ievTMP};
   %error('myprgm:outdate','this chose is no longer maintained')
  end
 end % i=1:numel(eigvecTMP)
else %wenn man den Vektor als vektor lässt, was eigentlich bullshit ist.
%  eigvec=eigvecTMP;
 for ievTMP=1:numel(eigvecTMP)
  eigvec{ievTMP}(:,:,1,:)=eigvecTMP{ievTMP}(:,:,:);
 end
end
EVsize=size(eigvec{1});
if model.numofeigs<=1
 s2=EVsize(2:end);
 sEig=model.numofeigs;
else 
 s2=EVsize(2:end-1);
 sEig=EVsize(end);
end
if realistic
 EVdofs=prod(EVsize(2:3));
else
 %assert(EVsize(1)==EVsize(2),'maybe should be symmetiric?')
 EVdofs=EVsize(2);
end
if forcedeig>EVsize(end)
 error('MyProgam:Input','forcedeig larger than num of eigenvalues')
end

%epsilon = model.lambda(2) - model.lambda(1);
%assert(epsilon-min(model.lambda(2:end) - model.lambda(1:end-1))<=+eps(4),'some steps are smaler than epsilon')
eigval = model.eigenvalues;

lambda0 = model.lambda0;
arclengths = model.arclengths;

fOld = length(eigval);
fMax =length(lambda0);
if fMax<fOld
 if fMax<52
  assert(fMax>=fOld,'set lambda0=5epsil or lager')
 end
 warning('MyPrgm:UnknownError','dont know what that means try forcerun=1')
 fOld=fMax;%try to fix assert
end
f=fOld;% f=min(fOld,find(lambda0==lastLambda));
eigposition = zeros(length(lambda0),1);

%% Computed Quantities:
V    = NaN(EVdofs,f); %V(:,k) = vs;
A    = NaN(EVdofs,f); %A(:,k) = as;
S    = NaN(f,1);            %S(k) = speeds;
A0   = NaN(f,1);            %A0(k) = a0s;
At   = NaN(f,1);            %At(k) = ats;
An   = NaN(f,1);            %An(k) = ans_;
RHO  = NaN(f,1);            %RHO(k) = rhos;
RHO2  = NaN(f,1);           %RHO2(k) = rhos2;
%RHO3  = zeros(f,1);           %RHO3(k) = rhos3;
%RHO4  = zeros(f,1);           %RHO4(k) = rhos4;
TAU  = NaN(f,1);            %TAU(k) = taus;
%TAU2  = zeros(f,1);           %TAU2(k) = tau2s;
OC0  = NaN(f,1);            %OC0(k) = ortConds0;
OC1  = NaN(f,1);            %OC1(k) = ortConds1;
OC2  = NaN(f,1);            %OC2(k) = ortConds2;
OC3  = NaN(f,1);            %OC3(k) = ortConds3;
OC4  = NaN(f,1);
OC5  = NaN(f,1);
OC6  = NaN(f,1);
OC7  = NaN(f,1);
OC8  = NaN(f,1);
OCeig  = NaN(f,1);          %OCeig(k) = 1;
LAM  = NaN(f,1);            %LAM(k) = lams;
EWd1l= NaN(f,1);
EWd2l= NaN(f,1);
oneEWd2l=NaN(f,1);
R    = NaN(EVdofs,f); %R(:,k) = r0s;
EVal= NaN(1,f);
r1ri = NaN(f,1);
% rsri = NaN(f,1);
drddr = NaN(f,1);
cosmu = NaN(f,1);
sinpsi = NaN(f,1);
coplanar = NaN(f,1);
%CosPhi = NaN(f,1);
X1 = NaN(f,1);X2 = NaN(f,1);X3 = NaN(f,1);X4 = NaN(f,1);
HYPO=NaN(f,1);
HypoB2110Zaeler=NaN(f,1);
CosGamma=NaN(f,1);
Normd3rds3=NaN(f,1);
SinGamma=NaN(f,1);
RXB=NaN(f,1);
DrhopDs=NaN(f,1);
Rconst=NaN(f,1);
EBENE=NaN(f,1);
PHIR=NaN(f,1);
%DetKt=NaN(f,1);
cosPhiMangA=NaN(f,1);
cosPhiMangB=NaN(f,1);
ZaelerB=NaN(f,1);
NormB1=NaN(f,1);
NormB2=NaN(f,1);
NormR=NaN(f,1);



% help quantities:
POS = NaN(f,1);             %POS(k) = poslam;
%%

lambda = lambda0(1:f);

[kl,lami] = InterpolateJK(model.fullEV,model.fulllambda);%[kl,lami] = findStabilityLimit(model.fullEV,model.fulllambda);
% if lami<0
%  kstability = kl;% kl - 1;
% else
%  kstability = kl;
% end
% if kstability == 0
%  kstability = 1;
% end
if kl==0
 kl=[];
end
if isempty(kl)
 warning('MyProgram:noStabilityLimit','assuming Stability-limit with last step')
 kl=numel(lambda0);
end
% stability_limit = lambda0(kl) + lami;
res.stability_limit = [kl,lami];
SchrittBefor=max(min(floor(res.stability_limit(1)),f),1);
EVatStability=eigvec{SchrittBefor}(evmiddle,:,:,1);
switch sortType
 case 'none'
  %disp('No sorting of eigenvectors');
 case 'forwards'
  %k = 1;
  poslam = 1;%; lami = 0;
 case 'backwards'
  [~,~,poslam] = findStabilityLimit(eigval,lambda0);
 otherwise
   poslam=forcedeig;
end
% kl=min([round(kl),numel(lambda0),numel(eigvec)]);
if ~strcmpi(sortType, 'none')
 %r0try = eigvec{k}(5,:,poslam);  %r0try = reshape(r0try,length(r0try),1);
 r1 = eigvec{1}(evmiddle,:,:,poslam);  r1 = reshape(r1,numel(r1),1);
%  rs = eigvec{kl}(evmiddle,:,:,poslam);
%  rs = reshape(rs,numel(rs),1);
else
 %r0try = zeros(size(eigvec{1},2),1);
 if ~isempty(forcedeig)
  r1 = eigvec{1}(evmiddle,:,:,forcedeig);
  r1 = reshape(r1,numel(r1),size(forcedeig,1));
%   rs = eigvec{kl}(evmiddle,:,:,forcedeig);
%   rs = reshape(rs,numel(rs),size(forcedeig,1));
 else
  r1 = eigvec{1}(evmiddle,:,:,1);  r1 = reshape(r1,numel(r1),1);
%   rs = eigvec{kl}(evmiddle,:,:,1);  rs = reshape(rs,numel(rs),1);
 end
end
r1=r1/norm(r1);
% rs=rs/norm(rs);

rho2atl0=NaN;
if ~isempty(forcedeig)
%  if numel(forcedeig)>0
  is0 = forcedeig;
%  else
%   is0=NaN;
%  end
else
 is0=NaN;
end

assert(numel(evmiddle)==1,'programming mistake');
s3alt=numel(is0);
tatl0=NaN([s2,s3alt]);
r0atl0=NaN([s2,s3alt]);

%rstabil=0.99992;%TL_arch3D-B31H-20-loadfac-1-eps0.01-KNL2-1.mat
%rstabil=0.99999981;%TL_arch3D-B31H-10-loadfac-1-eps0.01-KNL2-1.mat
%rstabil=0.9999999998;%TL_arch3D-B33H-10-loadfac-1-eps0.01-I-1.mat
%rstabil=0.999999999989;%TL_arch3D-B33H-10-loadfac-1-eps0.01-I-1.mat
%rstabil=0.999999999995;%TL_arch3D-B33H-10-loadfac-1-eps0.01-I-1.mat
%rstabil=0.9999999984;%TL_arch3D-B32-20-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
%rstabil=0.999999997;%TL_arch3D-B32H-20-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
%rstabil=0.9999999960;%TL_arch3D-B31H-10-loadfac-1-eps0.01-KNL2-1.mat (strengstens)
 if ~exist('main','var')
  main.undefined=true;
 end
 if ~isfield(main,'rstabil')
   main.rstabil=0.9999;
 end
rstabil=main.rstabil;

if strcmp(main.rho,'KtR1') || strcmp(main.Normierung,'rCT_K0_r') || strcmp(main.Normierung,'A0R1') || strcmp(main.whichEV,'sqrtK0_r') || strcmp(main.Normierung,'rK0r') || strcmp(main.Normierung,'k0_11')
 if sum(strcmp(fieldnames(model), 'stiffnessMatrices')) == 0
  error('MyPrgm:Missing','model.stiffnessMatrices missing, most likely main.whichEV changed, try modelprops.forcerun=1')
 end
 Kt0_0=model.stiffnessMatrices{1,1};
 if strcmp(main.whichEV,'sqrtK0_r')
  Kt0ii=full(diag(Kt0_0));
  dofs=uint16(numel(Kt0ii));
  Freiheitsgrade=7*dofs/10;% 7/10 to ignore hybrid dofs
  if any(Kt0ii(1:Freiheitsgrade)<=0)
   warning('MyPrgm:MangWrong','according to Mang this must be positive');
  end
  neg=find(Kt0ii<=0,1);
  if isempty(neg)
   removeDOF=inf;
  else
   removeDOF=max(neg,Freiheitsgrade+1);
   assert(removeDOF<uint16(inf),'overflow error')
  end
  if removeDOF<dofs
   Kt0ii(removeDOF:dofs)=0;
  end
  sqrtKt0ii=sqrt(Kt0ii); %only for main.whichEV='sqrtK0_r'
 end
end

if f<3
 warning('MyPrgm:strange','less than 3 steps')
end
for i = 1:f %f = length(eigval)

 eigposition(i) = is0(1);
 
 isi=1;%for isi = 1:1 %length(is)
 if ~any(isnan(r0atl0(:))) && strcmp(sortType,'forwardJK')
  for k=1:sEig
   rj=eigvec{i}(evmiddle,:,:,k);
   rj = reshape(rj,numel(rj),1);
   rj=rj/norm(rj);
   if norm(r0atl0-rj)<1 || norm(r0atl0+rj)<1 || norm(dot(r0atl0,rj))>0.33
    if is0(isi)~=k
     %warning('MyProgram:OderChange','The oder of %d eV %d changed to %d',forcedeig,is0(isi),k)
     is0(isi)=k;
    end
   end
  end
 end
 NrEw=is0(isi);
 if  strcmp(main.whichEV,'k0_11') || strcmp(main.whichEV,'k11') || strcmp(main.whichEV,'2023-12') || strcmp(main.whichEV,'2023_12half')
  r02 = model.eigvec2023{i}(:,1,NrEw);% r~(lambda-2*dlambda)
  r01 = model.eigvec2023{i}(:,2,NrEw);% r~ of previous loadstep
  rm  = model.eigvec2023{i}(:,3,NrEw);% r~ of current loadstep % DoFs x dlambda x NumberEigenValues
  r11 = model.eigvec2023{i}(:,4,NrEw);% r~ of next loadstep
  r12 = model.eigvec2023{i}(:,5,NrEw);% r~ of second next loadstep
 else
  r02 = eigvec{i}(evmiddle-2,:,:,NrEw);
  r01 = eigvec{i}(evmiddle-1,:,:,NrEw);
  rm = eigvec{i}(evmiddle,:,:,NrEw); % dlambda x Nodes x DoF-per-Node X NumberEigenvalues
  r11 = eigvec{i}(evmiddle+1,:,:,NrEw);
  r12 = eigvec{i}(evmiddle+2,:,:,NrEw);
 end
 r02 = reshape(r02,numel(r02),1);
 r02 = r02/norm(r02);
 r01 = reshape(r01,numel(r01),1);
 r01 = r01/norm(r01);
 if r02'*r01<0;  r01 = -r01;  end
 rm = reshape(rm,numel(rm),1);
 if strcmp(main.Normierung,'R1')
  if any(isnan(rm))
   warning('MyPrgm:rmNan:Sort','rm has NaNs')
   warning('off','MyPrgm:rmNan:Sort')
  end
  if i==3
   warning('on','MyPrgm:rmNan:Sort')
  end
  rm = rm/norm(rm);
 end
 eigvalSort=eigvalTMP{i}(evmiddle,NrEw);
 if r01'*rm<0;  rm = -rm;  end
 
 r11 = reshape(r11,numel(r11),1);
 r11 = r11/norm(r11);
 if rm'*r11<0;  r11 = -r11;  end
 r12 = reshape(r12,numel(r12),1);
 r12 = r12/norm(r12);
 if r11'*r12<0;  r12 = -r12;  end
 %if strcmp(main.Normierung,'k11') || strcmp(main.Normierung,'k0_11') ; eigvecH2i=model.eigvecH2{i}; end
 if strcmp(main.whichEV,'bungle_rKr') || strcmp(main.Normierung,'k11') || strcmp(main.whichEV,'k11')
  %rkDRH=squeeze(model.eigvecDRH{i}(1,:,:,NrEw));
  %rlDRH=squeeze(model.eigvecDRH{i}(2,:,:,NrEw));
  %rmDRH=squeeze(model.eigvecDRH{i}(3,:,:,NrEw)); % DoFpNode x Nodes
  %rnDRH=squeeze(model.eigvecDRH{i}(4,:,:,NrEw));
  %roDRH=squeeze(model.eigvecDRH{i}(5,:,:,NrEw));
  %Kt02=StiffMtxs{max(i-2,1)};
  Kt01=StiffMtxs{max(i-1,1)};
  KT  =StiffMtxs{i};
  Kt11=StiffMtxs{i+1};
%   if numel(StiffMtxs)<i+2
%    Kt12=NaN*Kt11;
%   else
%    Kt12=StiffMtxs{i+2};
%   end
 end
 
 
 if strcmp(main.whichEV,'bungle_rK0r') || strcmp(main.whichEV,'Disp_rK0r') || strcmp(main.Normierung,'rCT_K0_r') || strcmp(main.Normierung,'rK0r')
  Nenner02=sqrt(r02'*Kt0_0*r02);
  Nenner01=sqrt(r01'*Kt0_0*r01);
  Nenner0=sqrt(rm'*Kt0_0*rm);
  Nenner11=sqrt(r11'*Kt0_0*r11);
  Nenner12=sqrt(r12'*Kt0_0*r12);
 elseif strcmp(main.whichEV,'bungle_rKr')
  assert(0,'not implemented')
  Nenner01=sqrt(r01'*Kt01*r01);
  Nenner0=sqrt(rm'*KT*rm);
  Nenner11=sqrt(r11'*Kt11*r11);
 elseif strcmp(main.whichEV,'sqrtK_r') && strcmp(main.Normierung,'R1')
  iV=[i-2, i-1, i, i+1, i+2];
  iV(iV<1)=1;
  iV(iV>f)=f;
  r02=sqrt(full(diag(model.stiffnessMatrices{iV(1),1}))).*r02;
  r01=sqrt(full(diag(model.stiffnessMatrices{iV(2),1}))).*r01;
  rm=sqrt(full(diag(model.stiffnessMatrices{iV(3),1}))).*rm;
  r11=sqrt(full(diag(model.stiffnessMatrices{iV(4),1}))).*r11;
  r12=sqrt(full(diag(model.stiffnessMatrices{iV(5),1}))).*r12;
  Nenner02=norm(r02);
  Nenner01=norm(r01);
  Nenner0=norm(rm);
  Nenner11=norm(r11);
  Nenner12=norm(r12);
 elseif strcmp(main.whichEV,'sqrtK0_r') && strcmp(main.Normierung,'R1')
  r02=sqrtKt0ii.*r02;
  r01=sqrtKt0ii.*r01;
  rm=sqrtKt0ii.*rm;
  r11=sqrtKt0ii.*r11;
  r12=sqrtKt0ii.*r12;
  Nenner02=norm(r02);
  Nenner01=norm(r01);
  Nenner0=norm(rm);
  Nenner11=norm(r11);
  Nenner12=norm(r12);
 elseif strcmp(main.Normierung,'R1') ||  strcmp(main.Normierung,'k11') || strcmp(main.Normierung,'k0_11')
  if strcmp(main.Normierung,'k11')
   assert(strcmp(main.whichEV,'k11'),'this feature is not implemented')
  elseif strcmp(main.Normierung,'k0_11')
   assert(strcmp(main.whichEV,'k0_11'),'whichEV=%s and Normierung=%s  is not implemented',main.whichEV,main.Normierung)
  else
   assert(~strcmp(main.whichEV,'k11'),'this feature is not implemented')
   assert(~strcmp(main.whichEV,'k0_11'),'whichEV=%s and Normierung=%s is not implemented',main.whichEV,main.Normierung)
  end
  Nenner02=norm(r02);
  Nenner01=norm(r01);
  Nenner0=norm(rm);
  Nenner11=norm(r11);
  Nenner12=norm(r12);
 elseif strcmp(main.Normierung,'A0R1')
  Nenner02=norm(Kt0_0*r02);
  Nenner01=norm(Kt0_0*r01);
  Nenner0 =norm(Kt0_0*rm );
  Nenner11=norm(Kt0_0*r11);
  Nenner12=norm(Kt0_0*r12);
 elseif strcmp(main.Normierung,'skip')
  warning('MyPrgm:sortEigenV:Inconsistent:Input','Normierung is set to skip, but whichEV is not skip')
  warning('off','MyPrgm:sortEigenV:Inconsistent:Input')
  Nenner02=NaN;
  Nenner01=NaN;
  Nenner0 =NaN;
  Nenner11=NaN;
  Nenner12=NaN;
 else
  error('MyPrgm:NoTested','not tested/Implemented')
 end
 assert(numel(Nenner02)==1,'Matrix dimensions must agree.')
 r02 = r02/Nenner02;
 r01 = r01/Nenner01;
 rm = rm/Nenner0;
 r11 = r11/Nenner11;
 r12 = r12/Nenner12;


 if strcmp(main.rho,'A0R1')
  Mr02=Kt0_0*r02;
  Mr01=Kt0_0*r01;
  Mrm=Kt0_0*rm;
  Mr11=Kt0_0*r11;
  Mr12=Kt0_0*r12;
  if strcmp(main.Normierung,'A0R1')
   RS = [Mr02/norm(Mr02), Mr01/norm(Mr01),Mrm/norm(Mrm),Mr11/norm(Mr11),Mr12/norm(Mr12)];
  elseif  strcmp(main.Normierung,'R1')
   RS = [Mr02, Mr01,Mrm,Mr11,Mr12];
  else
   error('MyPrgm:NotTested','not tested/implemted')
  end
 elseif strcmp(main.rho,'R1')
  RS = [r02, r01, rm, r11, r12];
 elseif strcmp(main.rho,'skip')
  RS = [r02, r01, rm, r11, r12]*NaN;
 else
  error('MyPrgm:NotTested','not tested/implemted')
 end
 %disp(strcat('i=',num2str(i),': isnan(RS)=',num2str(all(isnan(RS(:))))))
 if ~all(isnan(RS(:)))
  dksi = arclengths{i};
  %            dksi = ones(size(dksi,1),size(dksi,2))*(lambda(2)-lambda(1));
  
  if any(isnan(r0atl0))
   %disp(['L508_forcedeig=',num2str(forcedeig),'; is0(isi)=',num2str(is0(isi)),'; i=',num2str(i)])
   if any(~isnan(r02)) && lambda0(i)~=0  && norm(dot(r02,r01))>rstabil
    r0atl0=r02;
   elseif any(~isnan(r01)) && lambda0(i)~=0  && norm(dot(r01,rm))>rstabil
    r0atl0=r01;
   elseif any(~isnan(rm))  && norm(dot(rm,r11))>rstabil
    if i<3%3 might be replaced by 2
     warning('MyProgram:strange','r0 at lambda0 must be NaN')
    end
    r0atl0=rm;
   elseif any(~isnan(r11))
    if norm(dot(r11,r12))>rstabil
     r0atl0=r11;
     %    elseif any(~isnan(r12)) && dot(r12,r13)>rstabil
     %     r0atl0=r12;
     %    elseif any(~isnan(r13)) && dot(r12,r13)>rstabil
     %     r0atl0=r13;
    else
     %warning('MyProgram:precission','r at beginning is imprecise')
    end
   else
    if i>2
     warning('MyProgram:strange','r11 at lambda0 schould not be NaN')
    end
   end
   %disp(r0atl0)
  end

  [v,a,speed,accel,at,an,rho,tau,ortCond1,rho2,drddr_,cosmu_,sinpsi_,ortCond2,ortCond3,ortCond4, t, ~, B, ~,...
   ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB,drhopds,rconst,Ebene,phiR,ortCond8,d2rds2] ...
   = getQuantities(RS,dksi,epsilon,r0atl0,rho2atl0,tatl0,main,EVatStability); % %DT=epsilon
  if isnan(rho2atl0) && i>1
   if abs(RHO2(i-1)-rho2)<.1
    rho2atl0=rho2;
   end
  end
  if any(isnan(tatl0(:))) && i>2
   tatl0=t;
  end
  NormB1(i)=norm(d2rds2);
  NormR(i)=norm(rm);
  if strcmp(main.whichEV,'bungle_rKr')
   ZaelerB(i)=(d2rds2'*Kt0_0*rm);
   NormB2(i)=norm(Kt0_0*rm);
  else
   ZaelerB(i)=NaN;
   NormB2(i)=NaN;
  end
  if any(isnan(t)) || ~strcmp(main.whichEV,'bungle_rKr')
   cosPhiMangA(i)=NaN;
   cosPhiMangB(i)=NaN;
  else
   cosPhiMangA(i)=-(t'*Kt0_0*t)/(norm(d2rds2)*norm(Kt0_0*rm));
   cosPhiMangB(i)=(d2rds2'*Kt0_0*rm)/(norm(d2rds2)*norm(Kt0_0*rm));
  end

  
  
  coplanar_ = B'*a/norm(a);
  %            %if rho>rhotry
  %rhotry = rho;
  %indi = 1;
  
  V(:,i) = v;
  A(:,i) = a;
  S(i) = speed;
  A0(i) = accel;
  At(i) = at;
  An(i) = an;
  RHO(i) = rho;
  if max(abs(imag(RHO)))>0.5
   warning('MyProgram:Complex','rho is komplex');
   warning('off','MyProgram:Complex')
  end
  RHO2(i) = rho2;
  %i
  TAU(i) = tau;
  OC0(i) = abs((v'*v)+(rm'*a));
  OC1(i) = ortCond1;
  OC2(i) = ortCond2;
  OC3(i) = ortCond3;
  OC4(i) = ortCond4;
  OC5(i) = ortCond5;
  OC6(i) = ortCond6;
  OC7(i) = ortCond7;
  OC8(i) = ortCond8;
  OCeig(i) = abs(rm'*v)/norm(v);
  %lami = ;
  LAM(i) = eigval{i}(5,is0(isi));
  EWd1l(i)= ( eigval{i}(6,is0(isi)) - eigval{i}(4,is0(isi)) ) / (2*epsilon);
  oneEWd2l(i)=1./EWd2l(i);
  %        disp(indi);
  R(:,i) = rm;%indi*rm;
  EVal(:,i)=eigvalSort;
  POS(i) = is0(1);
  
  r1ri(i) = abs(r1(:)'*rm);
%   rsri(i) = abs(rs(:)'*rm);
  drddr(i) = drddr_;
  cosmu(i) = cosmu_;
  sinpsi(i) = sinpsi_;
  coplanar(i) = coplanar_;
  %CosPhi(i)=cosPhi;
  X1(i)=x1;
  X2(i)=x2;
  X3(i)=x3;
  X4(i)=x4;
  HYPO(i)=Hypo;
  HypoB2110Zaeler(i)=NaN;%2*transpose(rm)*model.KB1{i}*rm;
  CosGamma(i)=cosGamma;
  Normd3rds3(i)=normd3rds3;
  SinGamma(i)=singamma;
  RXB(i)=RxB;
  DrhopDs(i)=drhopds;
  Rconst(i)=rconst;
  EBENE(i)=Ebene;
  PHIR(i)=phiR;
  
  %        r0try = indi*r0;
  %            end
 else
  warning('MyPrgm:RSNan:sortEigenValuesAndGetQuantities','RS is NAN (warning only shown once)')
  warning('off','MyPrgm:RSNan:sortEigenValuesAndGetQuantities')
  if i==13
   warning('on','MyPrgm:RSNan:sortEigenValuesAndGetQuantities')
  end
 end %all(isnan(RS))
 EWd2l(i)= ( eigval{i}(4,is0(isi)) - 2*eigval{i}(5,is0(isi)) +eigval{i}(6,is0(isi)) ) / (epsilon^2);
 LAM(i) = eigval{i}(5,is0(isi));
 %end %isi=1:1
end %i=1:f schleife ueber alle lambda

if any(RHO>1) && strcmp(main.Normierung,'R1')
 warning('MyProgramm:lowPrecission:RhoOne','rho is larger than one: %f',max(RHO))
 warning('off','MyProgramm:lowPrecission:RhoOne')
 RHO(RHO>1)=NaN;
else
 differences2=max(abs(RHO-RHO2));
 if differences2>2.004
  warning('MyProgramm:lowPrecission','rho differs from rho2 by %f',differences2)
 end
end
if any(RHO2>1) && strcmp(main.Normierung,'R1')
 warning('MyProgramm:lowPrecission:Rho2One','RHO2 is larger than one: %f',max(RHO2))
 warning('off','MyProgramm:lowPrecission:Rho2One')
 RHO2(RHO2>1)=NaN;
end
if any(RHO2<0)
 if any(RHO2<-1)
  warning('MyProgramm:lowPrecission:Rho2minusOne','RHO2 is smaller than -1: %f',min(RHO2))
  RHO2(RHO2<-1)=NaN;
 else
  warning('MyProgramm:lowPrecission:Rho2Zero','RHO2 is smaler than zero: %f (%d)',min(RHO2),sum(RHO2<0))
 warning('off','MyProgramm:lowPrecission:Rho2Zero')
  %RHO2(RHO2<0)=NaN;
 end
end

%% Apply orthogonality conditions

OC1(isnan(OC1))=0;
OC2(isnan(OC2))=0;
OC3(isnan(OC3))=0;
OC4(isnan(OC4))=0;
OC5(isnan(OC5))=0;
OC6(isnan(OC6))=0;
OC7(isnan(OC7))=0;
OC8(isnan(OC8))=0;

Orth1 = abs(OC1)<=limit.OC1; %0.0004;
Orth2 = abs(OC2)<=1e-12;%0.000001;
Orth3 = abs(OC3)<=3e-11;%0.000001;
Orth4 = abs(OC4)<=3e-13;%0.000001;
Orth5 = abs(OC5)<=limit.OC5; %=4.7e-5;% ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1 epsil = 0.01 % 3.6e-5%0.00006;
Orth5a = OC5<=limit.OC5a; %=0.0000044;% pureBendingBeam-B32OSH-2-len-5-loadfac--1 epsil = 0.01  %0.000003
Orth6 = abs(OC6)<=limit.OC6;
Orth6a = (OC6)<=limit.OC6a;
%Orth6a=Orth6;
Orth7 = abs(OC7)<=limit.OC7;
Orth7a = (OC7)<=eps(8);
Orth8 = abs(OC8)<.00561; %just to filter one result


if ~all(Orth5)
 warning('MyProgam:Limits:OC5abs','|OC5| is large: %f (%d)',max(abs(OC5)),sum(~Orth5))
  warning('off','MyProgam:Limits:OC5abs')
end
if ~all(Orth5a)
 warning('MyProgam:Limits:OC5pos','OC5 is positive: %f (%d)',max((OC5)),sum(~Orth5a))
  warning('off','MyProgam:Limits:OC5pos')
end
if ~all(Orth6)
 warning('MyProgam:Limits:OC6abs','|OC6| is large: %f (%d)',max(abs(OC6)),sum(~Orth6))
  warning('off','MyProgam:Limits:OC6abs')
end
if ~all(Orth6a)
 warning('MyProgam:Limits:OC6pos','OC6 is positive: %f (%d)',max((OC6)),sum(~Orth6a))
  warning('off','MyProgam:Limits:OC6pos')
end
if ~all(Orth7)
 warning('MyProgam:Limits:OC7large','OC7 is large: %f (%d)',max(abs(OC7)),sum(~Orth7))
  warning('off','MyProgam:Limits:OC7large')
end
if ~all(Orth7a)
 warning('MyProgam:Limits:OC7pos','OC7 is positive: %f',max((OC7)))
  warning('off','MyProgam:Limits:OC7pos')
end
if ~all(Orth8)
 warning('MyProgam:Limits:OC8','OC8 is large: %f (%d)',max(abs(OC8)),sum(~Orth8))
  warning('off','MyProgam:Limits:OC8')
end

Stmp=S;
Stmp(isnan(S))=0;
Cond2 = Stmp<limit.C2;
if ~all(Cond2)
 warning('MyProgam:Limits:Speed:Large','speed is large: %f',max(S))
  warning('off','MyProgam:Limits:Speed:Large')
end
A0(isnan(A0))=0;
Cond3 = A0<limit.C3A0; %790 % 1.28;%to be valid for c-KNL2-B32OS
if ~all(Cond3)
 warning('MyProgam:Limits:totAcc:Large','total accerlation is large: %f (%d)',max(A0),sum(~Cond3))
  warning('off','MyProgam:Limits:totAcc:Large')
end
At(isnan(At))=0;
Cond4 = abs(At)<limit.C4At; % pureBendingBeam-KNL2-B32OSH %1.09;%to be valid for c-KNL2-B32OS
if ~all(Cond4)
 warning('MyProgam:Limits:TanAcc:Large','tangential accerlation is large: %f (%d)',max(abs(At)),sum(~Cond4))
  warning('off','MyProgam:Limits:TanAcc:Large')
end
diffs=[At(2:end)-At(1:end-1);0];
Cond4Atdiff= (abs(diffs)<limit.C4Atdiff); % ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1 epsil = 0.01 %%%0.0090 %0.00070%1e-05
if ~all(Cond4Atdiff)
 warning('MyProgam:Limits:TanAcc:Fast','tangential accerlation diffs are changing fast: %f (%d)',max(abs(diffs)),sum(~Cond4Atdiff))
  warning('off','MyProgam:Limits:TanAcc:Fast')
end
An(isnan(An))=0;
Cond5An = abs(An)<limit.C5; % pureBendingBeam-KNL2-B32OSH 1.05;%to be valid for c-KNL2-B32OS 10Elements
if ~all(Cond5An)
 warning('MyProgam:Limits:NormalAcc:large','normal accerlation is large: %f (%d)',max(An),sum(~Cond5An))
  warning('off','MyProgam:Limits:NormalAcc:large')
end
diffs=[An(2:end)-An(1:end-1);0];
Cond5Andiff= (diffs>limit.C5diff);%-0.000017
if ~all(Cond5Andiff)
 warning('MyProgam:Limits:Acceleration:fast','normal accerlation diffs are changing fast: %f (%d)',min(diffs),sum(~Cond5Andiff))
  warning('off','MyProgam:Limits:Acceleration:fast')
end
%TAU(1)=0;
%TAU(isnan(TAU))=0;
Cond6 = abs(TAU)<limit.C6Tau;
%Cond6(3)=true;
if ~all(Cond6(3:end))
 warning('MyProgam:Limits:Tau:Large','TAU is large: %f (%d)',max(TAU),sum(~Cond6(3:end)))
  warning('off','MyProgam:Limits:Tau:Large')
end

%RHO2(isnan(RHO2))=0;
if strcmp(main.whichEV,'bungle_rK0r')
 Cond7 =true(f,1);
else
 Cond7 = ~(abs(RHO2)>1);
end
if ~all(Cond7)
 warning('MyProgam:Limits','|RHO2| is large: %f',max(abs(RHO2)))
end
Cond8 = ~(RHO2<limit.C8minrho);
if ~all(Cond8)
 warning('MyProgam:Limits:Rho2Small','RHO2 is small: %f (%d)',min(RHO2),sum(~Cond8))
 warning('off','MyProgam:Limits:Rho2Small')
end
tch=TAU(2:end)./TAU(1:end-1);
Cond9tmp = (tch<limit.C9) & (tch>1./limit.C9);
Cond9= [true;Cond9tmp] & [Cond9tmp;true];
if ~all(Cond9)
 warning('MyProgam:Limits:Tau:Change','tau changes by a factor of: %f (%d)',max(max(tch),1/min(tch)),sum(~Cond9))
  warning('off','MyProgam:Limits:Tau:Change')
end
 
OrthNaN=logical(Orth1.*Orth2.*Orth3.*Orth4.*Orth5.*Orth5a.*Orth6.*Orth6a.*Orth7.*Orth7a.*Orth8.*Cond2.*Cond3.*Cond4.*Cond4Atdiff.*Cond5An.*Cond5Andiff.*Cond6.*Cond7.*Cond8.*Cond9);
if ~all(OrthNaN)
 disp(strcat(num2str(numel(OrthNaN)),' values'))
end

%Orth=logical(Orth1);
if model.check==true
 %if false value set to NaN: find(~OrthNaN)
 %keep find(OrthNaN)
 %OrthNaN=logical(Orth1.*Orth2.*Orth3.*Orth4.*Orth5.*Orth5a.*Orth6.*Orth6a.*Orth7.*Orth7a.*Cond2.*Cond3.*Cond4.*Cond4Atdiff.*Cond5An.*Cond5Andiff.*Cond6.*Cond7.*Cond8);
else
 OrthNaN = true(length(Orth1),1);
end
Orth = true(length(Orth1),1);

res.lambdaorg = lambda;   
lambda(~OrthNaN)=NaN;
res.lambda = lambda(Orth)';
ZaelerB(~OrthNaN)=NaN;
res.ZaelerB=ZaelerB(Orth);
NormB1(~OrthNaN)=NaN;
res.NormB1=NormB1(Orth);
NormB2(~OrthNaN)=NaN;
res.NormB2=NormB2(Orth);
NormR(~OrthNaN)=NaN;
res.NormR=NormR(Orth);
cosPhiMangA(~OrthNaN)=NaN;
res.cosPhiMangA=cosPhiMangA(Orth);
cosPhiMangB(~OrthNaN)=NaN;
res.cosPhiMangB=cosPhiMangB(Orth);
V = V(:,Orth);
res.V = V;
A = A(:,Orth);            res.A = A;
S(~OrthNaN)=NaN;
res.S=S(Orth);
A0(~OrthNaN)=NaN;
res.A0 = A0(Orth);
At(~OrthNaN)=NaN;
res.At = At(Orth);
An(~OrthNaN)=NaN;
res.An = An(Orth);
res.RHO = RHO(Orth);
RHO2(~OrthNaN)=NaN;
res.RHO2 = RHO2(Orth);
%CosPhi(~OrthNaN)=NaN;
%res.CosPhi = CosPhi(Orth);
TAU(~OrthNaN)=NaN;
res.TAU = TAU(Orth);
%       OC1 = OC1(Orth);
OCeig = OCeig(Orth);      res.OCeig = OCeig;
LAM = LAM(Orth);          res.LAM = LAM;
EWd1l(~OrthNaN)=NaN;
res.EWd1l=real(EWd1l(Orth));
EWd2l(~OrthNaN)=NaN;
res.EWd2l=real(EWd2l(Orth));
oneEWd2l(abs(oneEWd2l)>1e2)=NaN;%filter out 1/zero
oneEWd2l(~OrthNaN)=NaN;
res.oneEWd2l=real(oneEWd2l(Orth));
R = R(:,Orth);
res.R = R;
EVal = EVal(:,Orth);
res.EVal = EVal;
POS = POS(Orth);          res.POS = POS;
r1ri = r1ri(Orth);        res.r1ri = r1ri;
% rsri = rsri(Orth);        res.rsri = rsri;
drddr = drddr(Orth);      res.drddr = drddr;
res.cosmu = cosmu(Orth);
Orthcosmu=[false;(res.cosmu(2:end)-res.cosmu(1:end-1))>0.6];
res.cosmu(Orthcosmu)=NaN;
sinpsi = sinpsi(Orth);    res.sinpsi = sinpsi;
coplanar = coplanar(Orth); res.coplanar = coplanar;
OC6(~OrthNaN)=NaN;
res.OC6 = OC6(Orth);
OC7(~OrthNaN)=NaN;
res.OC7 = OC7(Orth);

res.OC0 = OC0; res.OC1 = OC1; res.OC2 = OC2; res.OC3 = OC3; res.OC4 = OC4; res.OC5 = OC5;
res.eigposition = eigposition;
% res.stability_limit = [kl,stability_limit];

res.X1=X1(Orth);res.X2=X2(Orth);res.X3=X3(Orth);res.X4=X4(Orth);

res.HYPO=HYPO(Orth);

%res.HypoM2110=EWd2l./model.rddotKtr;
%res.HypoB2110=HypoB2110Zaeler./model.rddotKtr;
res.CosGamma=CosGamma(Orth);
res.Normd3rds3=Normd3rds3(Orth);
res.SinGamma=SinGamma(Orth);
RXB(~OrthNaN)=NaN;
res.RXB=RXB(Orth);
res.DrhopDs=DrhopDs(Orth);
res.Rconst=Rconst(Orth);
res.EBENE=EBENE(Orth);
res.PHIR=PHIR(Orth);

%res.rhorho3D = getQuantitiesStern(res.CosPhi,res.RHO2,res.lambda);
end
