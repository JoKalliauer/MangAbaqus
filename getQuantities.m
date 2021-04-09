function [v,a,dsdksi,accel,d2sdksi2,an,rho,tau,ortCond1,rho2, drddr, cosmu, sinpsi, ortCond2,ortCond3,ortCond4, t, N, B, r0,...
 ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB,drhopds,rconst,Ebene,phiR] ...
 = getQuantities(RS,~,DT,r0atl0,~,tatl0,main) %DT=epsilon
%rm4 = RS(:,1); %dksi04 = dksi(1);
%rm3 = RS(:,2); %dksi03 = dksi(2); ds04 = sqrt((r03-r04)'*(r03-r04));
rm2 = RS(:,3); %dksi02 = dksi(3); ds03 = sqrt((r02-r03)'*(r02-r03));
rm1 = RS(:,4); %dksi01 = dksi(4); ds02 = sqrt((r01-r02)'*(r01-r02));
r0  = RS(:,5);
rp1 = RS(:,6); %dksi11 = dksi(5); %ds01 = sqrt((r0-r01)'*(r0-r01));
rp2 = RS(:,7); %dksi12 = dksi(6); %ds11 = sqrt((r11-r0)'*(r11-r0));
%rp3 = RS(:,8); %dksi13 = dksi(7); %ds12 = sqrt((r12-r11)'*(r12-r11));
%rp4 = RS(:,9);
%dksi14 = dksi(8);
%ds13 = sqrt((r13-r12)'*(r13-r12));
%ds14 = sqrt((r14-r13)'*(r14-r13));

%r0check = (r11 + 3*r01)/(2*(0.5*(ds01+ds11))^2 + 4);

%    if nargin<3
%        DT = 0.25*(dksi11 + dksi12 + dksi13 + dksi14);
%    end
ortCond8=min([dot(rm2,r0),dot(rm1,r0),dot(rp1,r0),dot(rp2,r0)]);
%assert(ortCond8>.99,'the r vectors are to different')
NaNcheck=any(isnan(RS(:,4:6)));%any(isnan(RS(:,2:7)));
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
 %da=da/3;
 %      v = (r11 - r01)/(2);
 %      a = (r01 - 2*r0 + r11);
 %      da = (-1/2*r02 + 1*r01 - 1*r11 + 1/2*r12);
 
 %t02 = r01 - r03; %t02 = t02/norm(t02);
 %t01 = r0 - r02; %t01 = t01/norm(t01);
 %t0 = v/norm(v);
 %t11 = r12 - r0; %t11 = t11/norm(t11);
 %t12 = r13 - r11; %t12 = t12/norm(t12);
 
 %n01 = (t0 - t02);  %n01 = n01/norm(n01);
 %n0 = (t11 - t01);
 %n0 = n0/norm(n0);
 %n11 = (t12 - t0);  %n11 = n11/norm(n11);
 
 %t = t0;
 %n = n0;


accel = norm(a);

dsdksi = norm(v); %sd(13)
diff=abs(sqrt(v'*v)-dsdksi);
limit=3e-14;
if ~isnan(diff) && diff>limit
 if dsdksi>124
  warning('MyProgram:Boundary','Speed is large, returing partially NaN')
  v=v*NaN;
  a=a*NaN;
  r0=NaN*r0;
  dsdksi=NaN;
 else
  assert(diff<=limit,'norm(v) not equal sqrt(v*v) by %d',diff)
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
 warning('MyProgramm:NaN','cosGamma=%f is larger than one',cosGamma)
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
 warning('MyProgramm:lowPrecission','rho differs from 1/kappa by %s',differences)
end
if rho<0
 rho=NaN;
end
% if rho>0.8
%  disp(rho);
% end

absaXabsv=norm(a)*norm(v);
cosgamma=dot(v,a)/absaXabsv;
singamma=sqrt(1-cosgamma^2);

%rdxrdd=singamma*absaXabsv;

%tau20201113=norm(da)*cosgamma/rdxrdd;




%rsss=1/dsdksi^3*(1)
%rddd=(rp2 - 2*rp1 + 2*rm1 - rm2)/(2*(dsdksi*DT)^3);
%absrddd=norm(rddd);
normd3rds3=norm(d3rds3);

%dN = (kappa*d3rds3 - dkappa*d2rds2)/kappa^2;

 rho2 = -N'*r0;
 if rho2<0
  if rho2<-.1
   rho2=NaN;
  else
   warning('MyProgram:Output','rho2 is slighly negativ')
  end
 end
 differences=abs(rho-rho2);
 if differences>7.7939e-04%2.2088e-06
  %warning('MyProgramm:lowPrecission','rho differs from rho2 by %s',differences)
  if differences>.1 && main.check==true
   rho=NaN;
   rho2=NaN;
   warning('MyProgramm:lowPrecission','rho differs from rho2 by %s (to NaN)',differences)
  elseif main.check==true
   warning('MyProgramm:lowPrecission','rho differs from rho2 by %s (kept)',differences)
  end
 elseif isnan(differences) && main.check==true
  warning('MyProgramm:lowPrecission','rho or rho2 are NaN')
  rho=NaN;
  rho2=NaN;
 end

%tau2 = B'*dN;
RxB=dot(r0,B);
%rhodot=-dsdksi*tau*RxB;%-\dot{s}\,\tau\,(\mathrm{r}\cdot\mathrm{b})
drhopds=-tau*RxB;
rconst=r0atl0'*r0;
phiR=acos(abs(rconst));
Ebene=norm(r0-(tatl0'*r0)*tatl0-(r0'*r0atl0)*r0atl0);

 x1=1-rho^2/(1-rho^2)*drhopds^2;
 x2=dot(d3rds3,r0);
 x3=norm(rho2*d3rds3+t);
 x4=norm(r0-(tatl0'*r0)*tatl0-(r0'*r0atl0)*r0atl0);

Hypo=rho2^2*sqrt(1+tau^2);
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

 %    ortCond2 = drds'*d2rds2/norm(d2rds2);
 %    ortCond3 = drds'*d3rds3/norm(d3rds3);
 %    ortCond4 = d2rds2'*d3rds3/norm(d2rds2)/norm(d3rds3);
 
 %      tau2 = -drhodksi/(speed*(1 - rho^2)^(0.5));
 %      tau2 = -drhods/(norm(v)*(1 - rho^2)^(0.5));
 %      tau2 = norm(dn + kappa*t);
 %      tau2 = -drhodksi/(speed*(1 - rho^2)^(0.5));
 %      tau2 = -drho2/(speed*(B'*r0));
 %       tau2 = -drhodksi/(speed*(B'*r0));

 
 
 
 
 %    ortCond1 = abs((v'*v)./(r0'*a));
 ortCond1 = abs((r0'*a) + (v'*v));
 
 drddr = v'*a/norm(v)/norm(a);
 
 %cosmu = (r0'*a)/norm(a);
 cosmu = -((r0'*(rm1 - 2*r0 + rp1))/norm(rm1 - 2*r0 + rp1));
 cospsi = ((rp1 - rm1)'*(rp1 - 2*r0 + rm1))/(norm(rp1 - rm1)*norm(rp1 - 2*r0 + rm1));
 %sinpsi = sqrt(1 - (v'*a/norm(v)/norm(a)).^2);
 sinpsi = sqrt(1 - cospsi^2);
 %    rho = abs(cosmu/sinpsi);
 %R = r0;
 %V = v;
 %A = a;
 
 %by JoKalliauer
 ortCond6=r0'*a+dsdksi^2;
 
 ortCond5=r0'*t;
 
 % r=r0;
 % rdl=(r11-r01)/(2*DT);
 % rd2l=(r01-2*r0+r11)/(DT)^2;
 % ortCond6a=transpose(r)*rd2l+transpose(rdl)*rdl;
 ortCond7=2*transpose(r0)*(rm1+rp1)-transpose(rm1)*rp1-3;

 
 %tripJK=dot(r0,cross(T,N));
 %rr0=r0atl0'*r0;
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