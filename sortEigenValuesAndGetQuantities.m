function res = sortEigenValuesAndGetQuantities(model,sortType,~,forcedeig,limit,~,main) % forcedeig=k3;
if nargin<1
 model = runEigenProblem();
 sortType = 'backwards';
 %         sortType = 'forwards';
 %         sortType = 'none';
 %plotfig = [1,2];
 forcedeig = [];
end

epsilon = model.lambda(2) - model.lambda(1);
if ~exist('limit','var')
 limit.new=true;
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
  limit.OC5a=2.7e-5;
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
  limit.OC7=max(abs(-0.1482e-6),5.1600e-06*epsilon); %-5.1600e-06*epsilon from 
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
if strcmp(main.whichEV,'bungle')
 eigvecTMP = model.eigenvectors;
 realistic=false;
 evmiddle=5;
elseif strcmp(main.whichEV,'Disp')
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
elseif strcmp(main.whichEV,'all')
 eigvecTMP = model.eigvecDRH;
 realistic=true;
 evmiddle=3;
 relDofs=1:model.dofpNode;
 relNodes=1:model.inDOF(2);
else
%  eigvecTMP = model.eigvecDRH;
%  realistic=true;
%  evmiddle=3;
%  warning('MyProgam:whichEV','main.whichEV not regogniced')
 %relDofs=...;
 %relNodes=...;
 error('MyProgam:whichEV','main.whichEV not regogniced')
end
EVtmpsize=size(eigvecTMP{1});
eigvec=cell(size(eigvecTMP));
DimsEVtmp=numel(EVtmpsize);
if model.numofeigs<=1.5
 DimsEVtmp=DimsEVtmp+1;
end
if DimsEVtmp>3
%  eigvec=cell(size(eigvecTMP));
 for i=1:numel(eigvecTMP)
  if ~isempty(eigvecTMP{i})
   if isempty(relNodes)
    warning('MyProgram:Input','no relevant data found')
    res=[];
    return
    %forcedeig=[];
    %error('MyProgram:Input','no relevant data found')
   end
   eigvec{i}=eigvecTMP{i}(:,relDofs,relNodes,:);
  else
   eigvec{i}=eigvecTMP{i};
  end
 end
else
%  eigvec=eigvecTMP;
 for i=1:numel(eigvecTMP)
  eigvec{i}(:,:,1,:)=eigvecTMP{i}(:,:,:);
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

fold = length(eigval);
assert(length(lambda0)>=fold,'set lambda0=5epsil or lager')
f=fold;% f=min(fold,find(lambda0==lastLambda));
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
OCeig  = NaN(f,1);          %OCeig(k) = 1;
LAM  = NaN(f,1);            %LAM(k) = lams;
EWd2l= NaN(f,1);
oneEWd2l=NaN(f,1);
R    = NaN(EVdofs,f); %R(:,k) = r0s;
r1ri = NaN(f,1);
rsri = NaN(f,1);
drddr = NaN(f,1);
cosmu = NaN(f,1);
sinpsi = NaN(f,1);
coplanar = NaN(f,1);
%CosPhi = NaN(f,1);
X1 = NaN(f,1);X2 = NaN(f,1);X3 = NaN(f,1);X4 = NaN(f,1);
HYPO=NaN(f,1);
CosGamma=NaN(f,1);
Normd3rds3=NaN(f,1);
SinGamma=NaN(f,1);
RXB=NaN(f,1);
DrhopDs=NaN(f,1);
Rconst=NaN(f,1);
EBENE=NaN(f,1);
PHIR=NaN(f,1);
%DetKt=NaN(f,1);



% help quantities:
POS = NaN(f,1);             %POS(k) = poslam;
%%

lambda = lambda0(1:f);

[kl,lami] = findStabilityLimit(model.fullEV,model.fulllambda);
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
kl=min([round(kl),numel(lambda0),numel(eigvec)]);
if ~strcmpi(sortType, 'none')
 %r0try = eigvec{k}(5,:,poslam);  %r0try = reshape(r0try,length(r0try),1);
 r1 = eigvec{1}(evmiddle,:,:,poslam);  r1 = reshape(r1,numel(r1),1);
 rs = eigvec{kl}(evmiddle,:,:,poslam);  rs = reshape(rs,numel(rs),1);
else
 %r0try = zeros(size(eigvec{1},2),1);
 if ~isempty(forcedeig)
  r1 = eigvec{1}(evmiddle,:,:,forcedeig);
  r1 = reshape(r1,numel(r1),size(forcedeig,1));
  rs = eigvec{kl}(evmiddle,:,:,forcedeig);
  rs = reshape(rs,numel(rs),size(forcedeig,1));
 else
  r1 = eigvec{1}(evmiddle,:,:,1);  r1 = reshape(r1,numel(r1),1);
  rs = eigvec{kl}(evmiddle,:,:,1);  rs = reshape(rs,numel(rs),1);
 end
end
r1=r1/norm(r1);
rs=rs/norm(rs);

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


for i = 1:f %f = length(eigval)
 %disp(i)
%  C0 = eigvec{i}(5,:,:); C0 = reshape(C0,size(C0,2),size(C0,3));
%  
%  if strcmpi(sortType, 'none')
%   is0 = 1;
%  else
%   [~,is0] = compareEigenmodes(C0,r0try);
%  end
 %        val = [is0;val];
%  if ~isempty(forcedeig)
%   is0 = forcedeig;
%  end
 
 eigposition(i) = is0(1);
 %        Lami = eigval{i}(5,is);
 %        if sum(Lami<0)<length(Lami)
 %           is(Lami<0) = [];
 %        end
 
 %rhotry = 0;
 isi=1;%for isi = 1:1 %length(is)
 if ~any(isnan(r0atl0(:))) && strcmp(sortType,'forwardJK')
  for k=1:sEig
   rj=eigvec{i}(evmiddle,:,:,k);
   rj = reshape(rj,numel(rj),1);
   rj=rj/norm(rj);
   %rj = reshape(rj,length(rj),1);
   if norm(r0atl0-rj)<1 || norm(r0atl0+rj)<1 || norm(dot(r0atl0,rj))>0.33
    if is0(isi)~=k
     %warning('MyProgram:OderChange','The oder of %d eV %d changed to %d',forcedeig,is0(isi),k)
     is0(isi)=k;
    end
   end
  end
 end
 r02 = eigvec{i}(evmiddle-2,:,:,is0(isi));
 r02 = reshape(r02,numel(r02),1);
 r02 = r02/norm(r02);
 r01 = eigvec{i}(evmiddle-1,:,:,is0(isi)); r01 = reshape(r01,numel(r01),1);
 r01 = r01/norm(r01);
 if r02'*r01<0;  r01 = -r01;  end
 r0 = eigvec{i}(evmiddle,:,:,is0(isi));
 r0 = reshape(r0,numel(r0),1);
 r0 = r0/norm(r0);
 if r01'*r0<0;  r0 = -r0;  end
 r11 = eigvec{i}(evmiddle+1,:,:,is0(isi)); r11 = reshape(r11,numel(r11),1);
 r11 = r11/norm(r11);
 if r0'*r11<0;  r11 = -r11;  end
 r12 = eigvec{i}(evmiddle+2,:,:,is0(isi)); r12 = reshape(r12,numel(r12),1);
 r12 = r12/norm(r12);
 if r11'*r12<0;  r12 = -r12;  end
 %r13 = eigvec{i}(8,:,is0(isi)); r13 = reshape(r13,length(r13),1);
 %if r12'*r13<0;  r13 = -r13;  end
 %r14 = eigvec{i}(9,:,is0(isi)); r14 = reshape(r14,length(r14),1);
 %if r13'*r14<0;  r14 = -r14;  end
 

 
 
%  if (lambda0(kstability)-2*epsilon)<=lambda0(i) && (lambda0(i)<lambda0(kstability))
%   r11 = NaN*r11; r12 = NaN*r12; r13 = NaN*r13; r14 = NaN*r14;
%  elseif (lambda0(kstability)<lambda0(i))&&(lambda0(i)<=(lambda0(kstability)+2*epsilon))
%   r01 = NaN*r01; r02 = NaN*r02; r03 = NaN*r03; r04 = NaN*r04;
%  elseif lambda0(i)<2*epsilon
%   r01 = NaN*r01; r02 = NaN*r02; r03 = NaN*r03; r04 = NaN*r04;
%  end
 RS = [r02, r01, r0, r11, r12];%RS = [r04, r03, r02, r01, r0, r11, r12, r13, r14];
 if ~all(isnan(RS(:)))
  dksi = arclengths{i};
  %            dksi = ones(size(dksi,1),size(dksi,2))*(lambda(2)-lambda(1));
  
 if any(isnan(r0atl0))
  if any(~isnan(r02)) && lambda0(i)~=0  && norm(dot(r02,r01))>rstabil
   r0atl0=r02;
  elseif any(~isnan(r01)) && lambda0(i)~=0  && norm(dot(r01,r0))>rstabil
   r0atl0=r01;
  elseif any(~isnan(r0))  && norm(dot(r0,r11))>rstabil
   warning('MyProgram:strange','r0 at lambda0 must be NaN')
   r0atl0=r0;
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
   warning('MyProgram:strange','r11 at lambda0 schould not be NaN')
  end
  %disp(r0atl0)
 end

  [v,a,speed,accel,at,an,rho,tau,ortCond1,rho2,drddr_,cosmu_,sinpsi_,ortCond2,ortCond3,ortCond4, t, ~, B, ~,...
   ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB,drhopds,rconst,Ebene,phiR] ...
   = getQuantities(RS,dksi,epsilon,r0atl0,rho2atl0,tatl0,main); % %DT=epsilon
  %EWd1l= ( eigval{i}(4,is0(isi)) - eigval{i}(6,is0(isi)) ) / (2*epsilon);
  %x1=eigvec{i}(5,:,1)*(model.stiffnessMatrices{i,2}-EWd1l*model.stiffnessMatrices{1,1})*r0;
  if isnan(rho2atl0) && i>1
   if abs(RHO2(i-1)-rho2)<.1
    rho2atl0=rho2;
   end
  end
  if any(isnan(tatl0(:))) && i>2
   tatl0=t;
  end

  
  
  coplanar_ = B'*a/norm(a);
  %            %if rho>rhotry
  %rhotry = rho;
  indi = 1;
  
  V(:,i) = v;
  A(:,i) = a;
  S(i) = speed;
  A0(i) = accel;
  At(i) = at;
  An(i) = an;
  RHO(i) = rho;
  if max(abs(imag(RHO)))>0.5
   warning('MyProgram:Complex','rho is komplex');
  end
  RHO2(i) = rho2;
  TAU(i) = tau;
  OC0(i) = abs((v'*v)+(r0'*a));
  OC1(i) = ortCond1;
  OC2(i) = ortCond2;
  OC3(i) = ortCond3;
  OC4(i) = ortCond4;
  OC5(i) = ortCond5;
  OC6(i) = ortCond6;
  OC7(i) = ortCond7;
  OCeig(i) = abs(r0'*v)/norm(v);
  %lami = ;
  LAM(i) = eigval{i}(5,is0(isi));
  EWd2l(i)= ( eigval{i}(4,is0(isi)) - 2*eigval{i}(5,is0(isi)) +eigval{i}(6,is0(isi)) ) / (epsilon^2);
  oneEWd2l(i)=1./EWd2l(i);
  %        disp(indi);
  R(:,i) = indi*r0;
  POS(i) = is0(1);
  
  r1ri(i) = abs(r1(:)'*r0);
  rsri(i) = abs(rs(:)'*r0);
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
 end %all(isnan(RS))
 %end %isi=1:1
end %i=1:f schleife ueber alle lambda

% r1=RHO2.'.*cos(lambda.');
% r2=RHO2.'.*sin(lambda.');
% r3=sqrt(1-RHO2.'.^2);
% RS=[r1;r2;r3];
% res.rhoDach  = NaN(f,1);
% res.RXBDach =NaN(f,1);
% res.TauDach =NaN(f,1);
% res.SpeedDach=NaN(f,1);
% for i=5:f-4
%   %norm(RS(:,i))
%   [v,a,speed,accel,at,an,rho,tau,ortCond1,rho2,drddr_,cosmu_,sinpsi_,ortCond2,ortCond3,ortCond4, ~, ~, B, ~,...
%    ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB] ...
%    = getQuantities3D(RS(:,i-4:i+4),[],epsilon); %#ok<ASGLU>
%   res.rhoDach(i)=rho2;
%   res.RXBDach(i)=RxB;
%   res.TauDach(i) =tau;
%   res.SpeedDach(i)=speed;
% end

if any(RHO>1)
 warning('MyProgramm:lowPrecission','rho is larger than one: %f',max(RHO))
 RHO(RHO>1)=NaN;
else
 differences2=max(abs(RHO-RHO2));
 if differences2>2.004
  warning('MyProgramm:lowPrecission','rho differs from rho2 by %f',differences2)
 end
end
if any(RHO2>1)
 warning('MyProgramm:lowPrecission','RHO2 is larger than one: %f',max(RHO2))
 RHO2(RHO2>1)=NaN;
end
if any(RHO2<0)
 if any(RHO2<-1)
  warning('MyProgramm:lowPrecission','RHO2 is smaler than -1: %f',min(RHO2))
  RHO2(RHO2<-1)=NaN;
 else
  warning('MyProgramm:lowPrecission','RHO2 is smaler than zero: %f (%d)',min(RHO2),sum(RHO2<0))
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


if ~all(Orth5)
 warning('MyProgam:Limits','|OC5| is large: %f (%d)',max(abs(OC5)),sum(~Orth5))
end
if ~all(Orth5a)
 warning('MyProgam:Limits','OC5 is positive: %f (%d)',max((OC5)),sum(~Orth5a))
end
if ~all(Orth6)
 warning('MyProgam:Limits','|OC6| is large: %f (%d)',max(abs(OC6)),sum(~Orth6))
end
if ~all(Orth6a)
 warning('MyProgam:Limits','OC6 is positive: %f (%d)',max((OC6)),sum(~Orth6a))
end
if ~all(Orth7)
 warning('MyProgam:Limits','OC7 is large: %f (%d)',max(abs(OC7)),sum(~Orth7))
end
if ~all(Orth7a)
 warning('MyProgam:Limits','OC7 is positive: %f',max((OC7)))
end

S(isnan(S))=0;
Cond2 = S<limit.C2;
if ~all(Cond2)
 warning('MyProgam:Limits','speed is large: %f',max(S))
end
A0(isnan(A0))=0;
Cond3 = A0<limit.C3A0; %790 % 1.28;%to be valid for c-KNL2-B32OS
if ~all(Cond3)
 warning('MyProgam:Limits','total accerlation is large: %f (%d)',max(A0),sum(~Cond3))
end
At(isnan(At))=0;
Cond4 = abs(At)<limit.C4At; % pureBendingBeam-KNL2-B32OSH %1.09;%to be valid for c-KNL2-B32OS
if ~all(Cond4)
 warning('MyProgam:Limits','tangential accerlation is large: %f (%d)',max(abs(At)),sum(~Cond4))
end
diffs=[At(2:end)-At(1:end-1);0];
Cond4Atdiff= (abs(diffs)<limit.C4Atdiff); % ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1 epsil = 0.01 %%%0.0090 %0.00070%1e-05
if ~all(Cond4Atdiff)
 warning('MyProgam:Limits','tangential accerlation diffs are changing fast: %f (%d)',max(abs(diffs)),sum(~Cond4Atdiff))
end
An(isnan(An))=0;
Cond5An = abs(An)<limit.C5; % pureBendingBeam-KNL2-B32OSH 1.05;%to be valid for c-KNL2-B32OS 10Elements
if ~all(Cond5An)
 warning('MyProgam:Limits','normal accerlation is large: %f (%d)',max(An),sum(~Cond5An))
end
diffs=[An(2:end)-An(1:end-1);0];
Cond5Andiff= (diffs>limit.C5diff);%-0.000017
if ~all(Cond5Andiff)
 warning('MyProgam:Limits','normal accerlation diffs are changing fast: %f (%d)',min(diffs),sum(~Cond5Andiff))
end
%TAU(1)=0;
%TAU(isnan(TAU))=0;
Cond6 = abs(TAU)<limit.C6Tau;
%Cond6(3)=true;
if ~all(Cond6(3:end))
 warning('MyProgam:Limits','TAU is large: %f (%d)',max(TAU),sum(~Cond6(3:end)))
end

%RHO2(isnan(RHO2))=0;
Cond7 = ~(abs(RHO2)>1);
if ~all(Cond7)
 warning('MyProgam:Limits','|RHO2| is large: %f',max(abs(RHO2)))
end
Cond8 = ~(RHO2<limit.C8minrho);
if ~all(Cond8)
 warning('MyProgam:Limits','RHO2 is small: %f (%d)',min(RHO2),sum(~Cond8))
end
tch=TAU(2:end)./TAU(1:end-1);
Cond9tmp = (tch<limit.C9) & (tch>1./limit.C9);
Cond9= [true;Cond9tmp] & [Cond9tmp;true];
if ~all(Cond9)
 warning('MyProgam:Limits','tau changes by a factor of: %f (%d)',max(max(tch),1/min(tch)),sum(~Cond9))
end
 
OrthNaN=logical(Orth1.*Orth2.*Orth3.*Orth4.*Orth5.*Orth5a.*Orth6.*Orth6a.*Orth7.*Orth7a.*Cond2.*Cond3.*Cond4.*Cond4Atdiff.*Cond5An.*Cond5Andiff.*Cond6.*Cond7.*Cond8.*Cond9);
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
EWd2l(~OrthNaN)=NaN;
res.EWd2l=real(EWd2l(Orth));
oneEWd2l(abs(oneEWd2l)>1e2)=NaN;%filter out 1/zero
oneEWd2l(~OrthNaN)=NaN;
res.oneEWd2l=real(oneEWd2l(Orth));
R = R(:,Orth);
res.R = R;
POS = POS(Orth);          res.POS = POS;
r1ri = r1ri(Orth);        res.r1ri = r1ri;
rsri = rsri(Orth);        res.rsri = rsri;
drddr = drddr(Orth);      res.drddr = drddr;
res.cosmu = cosmu(Orth);
Orthcosmu=[false;(res.cosmu(2:end)-res.cosmu(1:end-1))>0.6];
res.cosmu(Orthcosmu)=NaN;
sinpsi = sinpsi(Orth);    res.sinpsi = sinpsi;
coplanar = coplanar(Orth); res.coplanar = coplanar;
%OC6(OC6<0)=0;% for debugging
OC6(~OrthNaN)=NaN;
res.OC6 = OC6(Orth);
%OC7(OC7<0)=0; %for debugging
OC7(~OrthNaN)=NaN;
res.OC7 = OC7(Orth);

res.OC0 = OC0; res.OC1 = OC1; res.OC2 = OC2; res.OC3 = OC3; res.OC4 = OC4; res.OC5 = OC5;
res.eigposition = eigposition;
% res.stability_limit = [kl,stability_limit];

res.X1=X1(Orth);res.X2=X2(Orth);res.X3=X3(Orth);res.X4=X4(Orth);

res.HYPO=HYPO(Orth);
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
