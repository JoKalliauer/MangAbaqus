function [v,a,dsdksi,accel,d2sdksi2,an,rho,tau,ortCond1,rho2, drddr, cosmu, sinpsi, ortCond2,ortCond3,ortCond4, t, N, B, r0,...
 ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB,drhopds,rconst,Ebene,phiR,ortCond8,d2rds2] ...
 = getQuantities(RS,~,DT,r0atl0,~,tatl0,main,EVatStability) %DT=epsilon
%% calulates tau, acceleration, rho, convergence-criterias
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)

%% Input
% RS  ... Eigenvectors
% dksi,~ ... removed
% DT ... epsilon-step-size, timeincrement
% r0atl0 ... rho at lambda=0 (not used in the current approach)
% ~ ...  removed
% tatl0 ... tau at lambda=0 (not used in the current approach)
% main ... structure with evaluation-settings from the start-file
% EVatStability ... Eigenvector at stabilitylimit  (not used in the current approach)


%% Output
% v ... Speed-vector of the eigenvector r
% a ... accelration-vector of the eigenvector r
% dsdksi ... ds/dlambda
% accel norm of  accerartion a
% d2sdksi2
% an .. normal accerlation
% rho ... eigenvector-curvature (outdated)
% tau .. torque of the vector-path
% ortCond1 .. orhogonalcodition check
% rho2 ... high precission eigenvector-curvature (better convergence)
% drddr 
% cosmu
% sinpsi
% ...

%% Recent Changes
%2023-02-17 JK: removed old comments and explained in/output
%2023-02-21 JK: removed commented out lines




%% Code

assert(size(RS,2)==5,'size of RS not 5 vectors');

rm2 = RS(:,1);%/norm(RS(:,1));
rm1 = RS(:,2);%/norm(RS(:,2)); 
r0  = RS(:,3);%/norm(RS(:,3));
rp1 = RS(:,4);%/norm(RS(:,4));
rp2 = RS(:,5);%/norm(RS(:,5)); 

ortCond8=min([dot(rm2,r0),dot(rm1,r0),dot(rp1,r0),dot(rp2,r0)]);

NaNcheck=any(isnan(RS(:,2:4)));%any(isnan(RS(:,2:7)));
if ~isfield(main,'rsame')
 main.rsame=0.8;
end
if ~isfield(main,'check')
 main.check=false %#ok<NOPRT>
end
if ortCond8<main.rsame && any(NaNcheck) && main.check==true
 whichNaNs=find(NaNcheck)+3;
 firstNaN=min(whichNaNs);
 lastNaN=max(whichNaNs);
 warning('MyProgram:changingR','there might be a change of order of eigenvektors, ortCond8=%f, RS is NaN:%d (%d:%d)',ortCond8,any(NaNcheck),firstNaN,lastNaN)
 if dot(rm2,r0)<main.rsame
  rm2=NaN*r0;
  %rm1=NaN*r0;
 end
 if dot(rm1,r0)<main.rsame
  rm1=NaN*r0;
 end
 if dot(rp1,r0)<main.rsame
  rp1=NaN*r0;
 end
 if dot(rp2,r0)<main.rsame
  rp2=NaN*r0;
  %rp1=NaN*r0;
 end
 r0=NaN*r0;
end


if ~isfield(main,'check')
 main.check=false %#ok<NOPRT>
end
if (~any(isnan(r0)) && ~any(isnan(rm1))) || main.check==false

 % Central difference
 v = (rp1 - rm1)/(2*DT); %(10)
 a = (rm1 - 2*r0 + rp1)/(DT)^2; %(11)
 da = 1/(DT)^3*(-1/2*rm2 + 1*rm1 - 1*rp1 + 1/2*rp2);%(12)

accel = norm(a);

dsdksi = norm(v); %sd; (13) in Report-January-2020.pdf; and (13) in Notes.pdf
diff=abs(sqrt(v'*v)-dsdksi);
limit=3e-14;
if ~isnan(diff) && diff>limit
 if dsdksi>124
  warning('MyProgram:Boundary:Speed','Speed is large, returing partially NaN')
  warning('off','MyProgram:Boundary:Speed')
  v=v*NaN;
  a=a*NaN;
  r0=NaN*r0;
  dsdksi=NaN;
 else
  assert(diff<=limit,'norm(v) not equal sqrt(v*v) by %d',diff)%does not work with single precission
 end
end
d2sdksi2 = (v'*a)/dsdksi;%sdd(14) or also at (tangential acceration)
d3sdksi3 = (((a'*a) + (v'*da))*dsdksi - (v'*a)*d2sdksi2)/(v'*v);%sddd(15)

t = v/dsdksi;%r'(16)
d2rds2 = 1/dsdksi^2*(a - d2sdksi2*t);% (17) r''=1/sd^2*(rdd-sdd*rd/sd)
d3rds3 = 1/dsdksi^3*(da - d3sdksi3*t - 3*dsdksi*d2sdksi2*d2rds2);%r'''(18)
%(da - d3sdksi3*drds - 3*dsdksi*d2sdksi2*d2rds2)=rddd-3*sdd*rdd/sd+(3*sdd^2/sd^2-sddd/sd)*rd
% - 3*dsdksi*d2sdksi2*d2rds2=-3*sdd*rdd/sd+3*sdd^2/sd^2*rd
% dsdksi*d2sdksi2*d2rds2=+sdd*rdd/sd-sdd^2/sd^2*rd
% sd*sdd*r''=+sdd*rdd/sd-sdd^2/sd^2*rd
% sd*sdd*(1/sd^2*(rdd-sdd*rd/sd))=+sdd*rdd/sd-(sdd^2/sd^2)*rd
% (rdd-sdd*rd/sd)=+rdd-(sdd/sd)*rd

%      speed = norm(v)/DT;
%      accel = norm(a)/DT^2;
%
%      dsdksi = sqrt(v'*v)/DT;
%      d2sdksi2 = (v'*a)/sqrt(v'*v)/DT^2;
%      d3sdksi3 = ((a'*a)/DT + v'*a)/(DT^2*sqrt(v'*v)) - ((v'*a)*(v'*a))/(DT^3*(v'*v)*sqrt(v'*v));
%
%      drds = v/sqrt(v'*v);
%      d2rds2 = a/(v'*v) - (v'*a)/(v'*v)^2*v;
%      d3rds3 = 1/dsdksi^3*(da - d3sdksi3*drds - 3*dsdksi*d2sdksi2*d2rds2);

kappa = norm(d2rds2); %Gl19
dkappa = (d3rds3'*d2rds2)/kappa;
%rho = 1/kappa;
%      rho = speed^3/sqrt(speed^2*accel^2 - (drds'*d2rds2)^2);

%drhodksi = -dsdksi*(d2rds2'*d3rds3)/(d2rds2'*d2rds2)^(3/2);
%drhods = -(d2rds2'*d3rds3)/(d2rds2'*d2rds2)^(3/2);
%      rho3 = rho; rho4 = rho;
tau = norm(1/kappa^2*(kappa*d3rds3 - dkappa*d2rds2) + t*kappa);
cosGamma=tau/norm(d3rds3);
if cosGamma>1
 warning('MyProgramm:cosGammaNaN','cosGamma=%f is larger than one',cosGamma)
 warning('off','MyProgramm:cosGammaNaN')
 cosGamma=NaN;
end





%at = d2sdksi2;
an = dsdksi^2 * kappa; %dsdksi^2/rho;

%T = drds;
N = d2rds2/kappa;
B = ((d3rds3*kappa - d2rds2*dkappa)/kappa^2 + t*kappa)/tau; %B = B/norm(B);

%rho = -N'*r0;
rho = norm(v)^3/sqrt((norm(v)*norm(a))^2 - (v'*a)^2);
differences=abs(rho-1/kappa);
if differences>2.2088e-06
 warning('MyProgramm:lowPrecissionKappa','rho differs from 1/kappa by %s',differences)
 warning('off','MyProgramm:lowPrecissionKappa')
end
if rho<0
 rho=NaN;
elseif rho>1.01
 warning('MyProgram:getQuanitites:OutputRHOlargeOver','rho is over 1.2')
 warning('off','MyProgram:getQuanitites:OutputRHOlargeOver')
 rho=NaN;
end
% if rho>0.8
%  disp(rho);
% end

absaXabsv=norm(a)*norm(v);

%cosgamma=dot(v,a)/absaXabsv;

%rdxrdd=singamma*absaXabsv;

%tau20201113=norm(da)*cosgamma/rdxrdd;




%rsss=1/dsdksi^3*(1)
%rddd=(rp2 - 2*rp1 + 2*rm1 - rm2)/(2*(dsdksi*DT)^3);
%absrddd=norm(rddd);
normd3rds3=norm(d3rds3);

%dN = (kappa*d3rds3 - dkappa*d2rds2)/kappa^2;

 if strcmp(main.Normierung,'rCT_K0_r')
  rho2=rho;
 else
  rho2 = -N'*r0;
 end
 if rho2<0
  if rho2<-.01
   if strcmp(main.Normierung,'R1')
    rho2=NaN;
   elseif strcmp(main.Normierung,'A0R1')
    rho2=NaN; 
   elseif abs(norm(r0)-1)<=eps(1)
    warning('MyPrgm:NotTested','not tested')
   else
    warning('MyPrgm:NotTested','not tested')
   end
  else
   warning('MyProgram:OutputRHONegative','rho2 is slighly negativ')
   warning('off','MyProgram:OutputRHONegative')
  end
 end
 differences=abs(rho-rho2);
 if differences>7.7939e-04%2.2088e-06
  %warning('MyProgramm:lowPrecission','rho differs from rho2 by %s',differences)
  if differences>.1 && main.check==true
   rho=NaN;
   rho2=NaN;
   warning('MyProgramm:lowPrecissionRHONaN','rho differs from rho2 by %s (to NaN)',differences)
   warning('off','MyProgramm:lowPrecissionRHONaN')
  elseif main.check==true
   warning('MyProgramm:lowPrecissionRHOKept','rho differs from rho2 by %s (kept)',differences)
   warning('off','MyProgramm:lowPrecissionRHOKept')
  end
 elseif isnan(differences) && main.check==true
  warning('MyProgramm:lowPrecissionRHO2','rho or rho2 are NaN')
  warning('off','MyProgramm:lowPrecissionRHO2')
  rho=NaN;
  rho2=NaN;
 end

%tau2 = B'*dN;
RxB=dot(r0,B);
%rhodot=-dsdksi*tau*RxB;%-\dot{s}\,\tau\,(\mathrm{r}\cdot\mathrm{b})
drhopds=-tau*RxB;
if exist('EVatStability','var') && strcmp(main.whichEV,'bungle')
 rconst=dot(EVatStability,r0);%(EVatStability*r0);
else
 rconst=r0atl0(:)'*r0;
end
phiR=acos(abs(rconst));
Ebene=norm(r0-(tatl0(:)'*r0)*tatl0(:)-(r0'*r0atl0(:))*r0atl0(:));

 %x1=dot(r0,d2rds2)+1;
 %x2=1/norm(d2rds2);
 %x3=norm(rho2*d3rds3+t);
 %x4=abs(dot(r0,d2rds2)+1);
 
 x1=abs(transpose(t)*d2rds2);
 x2=r0'*d2rds2;
 x3=norm(rho2*d3rds3+t);
 x4=abs(dot(r0,d2rds2)+1);

 Hypo=rho2*(1+tau);
 if imag(Hypo)~=0
  Hypo=NaN;
 end
 
 ortCond2 = t'*N;
 if ~isnan(ortCond2)
  if abs(ortCond2)>1e-12
   %assert(abs(ortCond2)<=1e-12,'ortCond2 not fullfiled by %d',abs(ortCond2))
  end
 end
 ortCond3 = t'*B;
 ortCond4 = N'*B;
 ortCond1 = abs((r0'*a) + (v'*v));
 
 drddr = dot(v,a)/(absaXabsv);
 singamma=sqrt(1-drddr^2);

 cosmu = -((r0'*(rm1 - 2*r0 + rp1))/norm(rm1 - 2*r0 + rp1));
 cospsi = ((rp1 - rm1)'*(rp1 - 2*r0 + rm1))/(norm(rp1 - rm1)*norm(rp1 - 2*r0 + rm1));
 sinpsi = sqrt(1 - cospsi^2);

 %by JoKalliauer
 ortCond6=r0'*a+dsdksi^2;
 
 ortCond5=r0'*t;
 ortCond7=2*transpose(r0)*(rm1+rp1)-transpose(rm1)*rp1-3;
 ortCond8=dot(r0,d2rds2)+1;
else
 warning('MyProgram:NaN','returning NaN')
 v=NaN*r0;
 a=NaN*r0;
 dsdksi=NaN;
 accel=NaN;
 d2sdksi2=NaN;
 an=NaN;
 rho=NaN;
 tau=NaN;
 ortCond2=NaN;
 ortCond3=NaN;
 ortCond4=NaN;
 t=NaN;
 d2rds2=NaN*r0;
 N=NaN;
 B=NaN*r0;
 Hypo=NaN;
 cosGamma=NaN;
 normd3rds3=NaN;
 singamma=NaN;
 x1=NaN;
 x2=NaN;
 x3=NaN;
 x4=NaN;
 RxB=NaN;
 drhopds=NaN;
 rconst=NaN;
 phiR=NaN;
 Ebene=NaN;
 
 rho2=NaN;
 ortCond1=NaN;
 drddr=NaN;
 cosmu=NaN;
 sinpsi=NaN;
 ortCond6=NaN;
 ortCond5=NaN;
 ortCond7=NaN;
 %tripJK=NaN;
end

% if exist('r0atl0','var') && exist('rho2atl0','var')
%  if ~isnan(rho2)
%   cosTheta=sqrt(1-rho2^2);
%   cosTheta0=sqrt(1-rho2atl0^2);
%   Theta=acos(cosTheta);
%   Theta0=acos(cosTheta0);
%  else
%   cosTheta=NaN;
%   cosTheta0=NaN;
%   Theta=NaN;
%   Theta0=NaN;
%  end
%  x1=abs(r0atl0.'*r0)*cos(Theta-Theta0);
%  x2=cosTheta0*cosTheta;
%  x3=(x1- x2 );
%  x4=(rho2atl0*rho2);
%  cosPhi= x3 /x4;
%  if abs(cosPhi)-1>1.8
%   %warning('MyProgram:Inpossible','cosPhi might be between -1 and 1')
%  end
% end
end