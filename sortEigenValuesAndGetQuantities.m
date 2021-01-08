function res = sortEigenValuesAndGetQuantities(model,sortType,~,forcedeig,limit,~)
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
 limit.OC6a=5e-12;
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
 limit.C6Tau=3.9041e7;
end


if ~isfield(limit,'C8minrho')
 limit.C8minrho=0;
end


if sum(strcmp(fieldnames(model), 'check')) == 0
 model.check=true;
end

eigvec = model.eigenvectors;
if forcedeig>size(eigvec{1},3)
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
V    = NaN(size(eigvec{1},2),f); %V(:,k) = vs;
A    = NaN(size(eigvec{1},2),f); %A(:,k) = as;
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
R    = NaN(size(eigvec{1},2),f); %R(:,k) = r0s;
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
% help quantities:
POS = NaN(f,1);             %POS(k) = poslam;
%%

lambda = lambda0(1:f);

[kl,lami,~] = findStabilityLimit(eigval,lambda0);
% if lami<0
%  kstability = kl;% kl - 1;
% else
%  kstability = kl;
% end
% if kstability == 0
%  kstability = 1;
% end

if kl==0
 warning('MyProgram:noStabilityLimit','assuming Stability-limit with last step')
 kl=numel(lambda0);
end
stability_limit = lambda0(kl) + lami;
switch sortType
 case 'none'
  %disp('No sorting of eigenvectors');
 case 'forwards'
  k = 1; poslam = 1;%; lami = 0;
 case 'backwards'
  [k,~,poslam] = findStabilityLimit(eigval,lambda0);
end

if ~strcmpi(sortType, 'none')
 r0try = eigvec{k}(5,:,poslam);  r0try = reshape(r0try,length(r0try),1);
 r1 = eigvec{1}(5,:,poslam);  r1 = reshape(r1,length(r1),1);
 rs = eigvec{kl}(5,:,poslam);  rs = reshape(rs,length(rs),1);
else
 r0try = zeros(size(eigvec{1},2),1);
 if ~isempty(forcedeig)
  r1 = eigvec{1}(5,:,forcedeig);
  r1 = reshape(r1,length(r1),size(r1,3));
  rs = eigvec{kl}(5,:,forcedeig);
  rs = reshape(rs,length(rs),size(rs,3));
 else
  r1 = eigvec{1}(5,:,1);  r1 = reshape(r1,length(r1),1);
  rs = eigvec{kl}(5,:,1);  rs = reshape(rs,length(rs),1);
 end
end

rho2atl0=NaN;
[s1,s2,s3]=size(eigvec{1}(5,:,:));
assert(s1==1,'programming mistake');
tatl0=NaN(s2,1);
for i = 1:f %f = length(eigval)
 %disp(i)
 C0 = eigvec{i}(5,:,:); C0 = reshape(C0,size(C0,2),size(C0,3));
 
 if strcmpi(sortType, 'none')
  is0 = 1;
 else
  [~,is0] = compareEigenmodes(C0,r0try);
 end
 %        val = [is0;val];
 if ~isempty(forcedeig)
  is0 = forcedeig;
 end
 
 eigposition(i) = is0(1);
 %        Lami = eigval{i}(5,is);
 %        if sum(Lami<0)<length(Lami)
 %           is(Lami<0) = [];
 %        end
 
 %rhotry = 0;
 isi=1;%for isi = 1:1 %length(is)
 r04 = eigvec{i}(1,:,is0(isi)); r04 = reshape(r04,length(r04),1);
 r03 = eigvec{i}(2,:,is0(isi)); r03 = reshape(r03,length(r03),1);
 if r04'*r03<0;  r03 = -r03;  end
 r02 = eigvec{i}(3,:,is0(isi)); r02 = reshape(r02,length(r02),1);
 if r03'*r02<0;  r02 = -r02;  end
 r01 = eigvec{i}(4,:,is0(isi)); r01 = reshape(r01,length(r01),1);
 if r02'*r01<0;  r01 = -r01;  end
 r0 = eigvec{i}(5,:,is0(isi)); r0 = reshape(r0,length(r0),1);
 if r01'*r0<0;  r0 = -r0;  end
 r11 = eigvec{i}(6,:,is0(isi)); r11 = reshape(r11,length(r11),1);
 if r0'*r11<0;  r11 = -r11;  end
 r12 = eigvec{i}(7,:,is0(isi)); r12 = reshape(r12,length(r12),1);
 if r11'*r12<0;  r12 = -r12;  end
 r13 = eigvec{i}(8,:,is0(isi)); r13 = reshape(r13,length(r13),1);
 if r12'*r13<0;  r13 = -r13;  end
 r14 = eigvec{i}(9,:,is0(isi)); r14 = reshape(r14,length(r14),1);
 if r13'*r14<0;  r14 = -r14;  end
 

 
 
%  if (lambda0(kstability)-2*epsilon)<=lambda0(i) && (lambda0(i)<lambda0(kstability))
%   r11 = NaN*r11; r12 = NaN*r12; r13 = NaN*r13; r14 = NaN*r14;
%  elseif (lambda0(kstability)<lambda0(i))&&(lambda0(i)<=(lambda0(kstability)+2*epsilon))
%   r01 = NaN*r01; r02 = NaN*r02; r03 = NaN*r03; r04 = NaN*r04;
%  elseif lambda0(i)<2*epsilon
%   r01 = NaN*r01; r02 = NaN*r02; r03 = NaN*r03; r04 = NaN*r04;
%  end
 RS = [r04, r03, r02, r01, r0, r11, r12, r13, r14];
 if ~all(isnan(RS(:)))
  dksi = arclengths{i};
  %            dksi = ones(size(dksi,1),size(dksi,2))*(lambda(2)-lambda(1));
  
 if ~exist('r0atl0','var')
   if any(~isnan(r04)) && lambda0(i)~=0
    r0atl0=r04;
   elseif any(~isnan(r03)) && lambda0(i)~=0
    r0atl0=r03;
   elseif any(~isnan(r02)) && lambda0(i)~=0
    r0atl0=r02;
   elseif any(~isnan(r01)) && lambda0(i)~=0
    r0atl0=r01;
   elseif any(~isnan(r0))
    warning('MyProgram:strange','r0 at lambda0 must be NaN')
    r0atl0=r0;
   elseif any(~isnan(r11))
    r0atl0=r11;
   elseif any(~isnan(r12))
    r0atl0=r12;
   elseif any(~isnan(r13))
    r0atl0=r13;
   else
    warning('MyProgram:strange','r11 at lambda0 schould not be NaN')
   end
 end

  [v,a,speed,accel,at,an,rho,tau,ortCond1,rho2,drddr_,cosmu_,sinpsi_,ortCond2,ortCond3,ortCond4, t, ~, B, ~,...
   ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB] ...
   = getQuantities(RS,dksi,epsilon,r0atl0,rho2atl0,tatl0); %DT=epsilon
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
  
  r1ri(i) = abs(r1'*r0);
  rsri(i) = abs(rs'*r0);
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
  warning('MyProgramm:lowPrecission','RHO2 is smaler than zero: %f',min(RHO2))
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
 warning('MyProgam:Limits','OC6 is large: %f (%d)',max(abs(OC6)),sum(~Orth6))
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
Cond6(3)=true;
if ~all(Cond6(3:end))
 warning('MyProgam:Limits','TAU is large: %f',max(TAU))
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
 
OrthNaN=logical(Orth1.*Orth2.*Orth3.*Orth4.*Orth5.*Orth5a.*Orth6.*Orth6a.*Orth7.*Orth7a.*Cond2.*Cond3.*Cond4.*Cond4Atdiff.*Cond5An.*Cond5Andiff.*Cond6.*Cond7.*Cond8);
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
res.stability_limit = [kl,stability_limit];

res.X1=X1(Orth);res.X2=X2(Orth);res.X3=X3(Orth);res.X4=X4(Orth);

res.HYPO=HYPO(Orth);
res.CosGamma=CosGamma(Orth);
res.Normd3rds3=Normd3rds3(Orth);
res.SinGamma=SinGamma(Orth);
res.RXB=RXB(Orth);

%res.rhorho3D = getQuantitiesStern(res.CosPhi,res.RHO2,res.lambda);
end
