function [v,a,dsdksi,accel,d2sdksi2,an,rho,tau,ortCond1,rho2, drddr, cosmu, sinpsi, ortCond2,ortCond3,ortCond4, t, N, B, r0,...
 ortCond5,ortCond6,ortCond7,Hypo,cosGamma,normd3rds3,singamma,x1,x2,x3,x4,RxB] ...
 = getQuantities3D(RS,~,DT,~,~) %DT=epsilon
%rm4 = RS(:,1); %dksi04 = dksi(1);
%rm3 = RS(:,2); %dksi03 = dksi(2); ds04 = sqrt((r03-r04)'*(r03-r04));
rm2 = RS(:,3); %dksi02 = dksi(3); ds03 = sqrt((r02-r03)'*(r02-r03));
rm1 = RS(:,4); %dksi01 = dksi(4); ds02 = sqrt((r01-r02)'*(r01-r02));
r0  = RS(:,5);
rp1 = RS(:,6); %dksi11 = dksi(5); %ds01 = sqrt((r0-r01)'*(r0-r01));
rp2 = RS(:,7); %dksi12 = dksi(6); %ds11 = sqrt((r11-r0)'*(r11-r0));
%rp3 = RS(:,8); %dksi13 = dksi(7); %ds12 = sqrt((r12-r11)'*(r12-r11));
%rp4 = RS(:,9);

 v = (rp1 - rm1)/(2*DT); %(10)
 a = (rm1 - 2*r0 + rp1)/(DT)^2; %(11)
 da = 1/(DT)^3*(-1/2*rm2 + 1*rm1 - 1*rp1 + 1/2*rp2);%(12)

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

kappa = norm(d2rds2); %Gl19
%dkappa = (d3rds3'*d2rds2)/kappa;

VxA=cross(v,a);
tau=dot(VxA,da)/dot(VxA,VxA);
%tau = norm(1/kappa^2*(kappa*d3rds3 - dkappa*d2rds2) + t*kappa);
cosGamma=tau/norm(d3rds3);
if cosGamma>1
 warning('MyProgramm:NaN','cosGamma=%f is larger than one',cosGamma)
 cosGamma=NaN;
end





%at = d2sdksi2;
an = dsdksi^2 * kappa; %dsdksi^2/rho;

%T = drds;
N = d2rds2/kappa;
%B = ((d3rds3*kappa - d2rds2*dkappa)/kappa^2 + t*kappa)/tau; 
absaXabsv=norm(a)*norm(v);
cosgamma=dot(v,a)/absaXabsv;
singamma=sqrt(1-cosgamma^2);
B=cross(v,a)/absaXabsv;
B = B/norm(B);

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




normd3rds3=norm(d3rds3);


RxB=dot(r0,B);

drhopds=-tau*RxB;
 x1=NaN;
 x2=1-rho^2/(1-rho^2)*drhopds^2;
 x3=RxB;%3*dot(a,v)+dot(r0,da);
 x4=-tau*RxB;

Hypo=rho*(1+RxB^2);
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
if ~any(isnan(r0))

 rho2 = -N'*r0;
 if rho2<0
  if rho2<-.1
   rho2=NaN;
  else
   warning('MyProgram:Output','rho2 is slighly negativ')
  end
 end
 
 
 
 

 ortCond1 = abs((r0'*a) + (v'*v));
 
 drddr = v'*a/norm(v)/norm(a);
 

 cosmu = -((r0'*(rm1 - 2*r0 + rp1))/norm(rm1 - 2*r0 + rp1));
 cospsi = ((rp1 - rm1)'*(rp1 - 2*r0 + rm1))/(norm(rp1 - rm1)*norm(rp1 - 2*r0 + rm1));

 sinpsi = sqrt(1 - cospsi^2);

 

 ortCond6=r0'*a+dsdksi^2;
 
 ortCond5=r0'*t;
 

 ortCond7=2*transpose(r0)*(rm1+rp1)-transpose(rm1)*rp1-3;
 

else
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

end