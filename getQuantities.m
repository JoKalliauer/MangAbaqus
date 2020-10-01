function [v,a,speed,accel,at,an,rho,tau,ortCond1,rho2, drddr, cosmu, sinpsi, ortCond2,ortCond3,ortCond4, T, N, B, R, V, A] ...
 = getQuantities(RS,dksi,DT) %DT=epsilon
   r04 = RS(:,1); %dksi04 = dksi(1); 
   r03 = RS(:,2); %dksi03 = dksi(2); ds04 = sqrt((r03-r04)'*(r03-r04));
   r02 = RS(:,3); %dksi02 = dksi(3); ds03 = sqrt((r02-r03)'*(r02-r03));
   r01 = RS(:,4); %dksi01 = dksi(4); ds02 = sqrt((r01-r02)'*(r01-r02));
   r0  = RS(:,5);
   r11 = RS(:,6); dksi11 = dksi(5); %ds01 = sqrt((r0-r01)'*(r0-r01));
   r12 = RS(:,7); dksi12 = dksi(6); %ds11 = sqrt((r11-r0)'*(r11-r0));
   r13 = RS(:,8); dksi13 = dksi(7); %ds12 = sqrt((r12-r11)'*(r12-r11));
   r14 = RS(:,9); dksi14 = dksi(8); %ds13 = sqrt((r13-r12)'*(r13-r12));
                                    %ds14 = sqrt((r14-r13)'*(r14-r13));
                                    
   %r0check = (r11 + 3*r01)/(2*(0.5*(ds01+ds11))^2 + 4);
   
   if nargin<3
       DT = 0.25*(dksi11 + dksi12 + dksi13 + dksi14);
   end
   if sum(abs(r01))==0 % Forward difference
%     % FIRST ORDER
%      v = 1/dksi11*(-r0 + r11);
%      a = 1/(0.5*(dksi11 + dksi12))^2*(r0 -2*r11 + r12);
%      da = 1/(1/3*(dksi11 + dksi12 + dksi13))^3*(-1*r0 + 3*r11 - 3*r12 +1*r13);
% 
%      t0 = v/norm(v);
%      t11 = r11 - r0; t11 = t11/norm(t11);
%      t12 = r12 - r11; t12 = t12/norm(t12);
%        
%      n0 = (t11 - t0);  
%      n0 = n0/norm(n0);
%      n11 = (t12 - t11);  n11 = n11/norm(n11);
%      
%      t = t0;
%      n = n0;
%      
%      dn = (n11 - n0)/(ds11);
%           
%      rho2 = -n'*r0;
%      rho11 = -n11'*r11;
%      drho2 = 1/(dksi11)*(rho11 - rho2);
% %      drho2 = 1/(ds11)*(rho11 - rho2);

     % SECOND ORDER:
     v = 1/(2*DT)*(-r12 + 4*r11 - 3*r0);
     a = 1/(DT)^2*(2*r0 - 5*r11 + 4*r12 - 1*r13);
     da = 1/DT^3*(-5/2*r0 + 9*r11 - 12*r12 + 7*r13 - 3/2*r14);
%      v = 1/2*(-r12 + 4*r11 - 3*r0);
%      a = (2*r0 - 5*r11 + 4*r12 - 1*r13);
%      da = (-5/2*r0 + 9*r11 - 12*r12 + 7*r13 - 3/2*r14);

     %t0 = v/norm(v);
     %t11 = r12 - r0; %t11 = t11/norm(t11);
     %t12 = r13 - r11; %t12 = t12/norm(t12);
     %t13 = r14 - r12; %t13 = t13/norm(t13);
       
     %n0 = (-t12 + 4*t11 - 3*t0);  
     %n0 = n0/norm(n0);
     %n11 = (t12 - t0);  %n11 = n11/norm(n11);
     %n12 = (t13 - t11);  %n12 = n12/norm(n12);
     
     %t = t0;
     %n = n0;
          
   elseif sum(abs(r11))==0 % Backward difference
     % SECOND ORDER:
     v = 1/(2*DT)*(r02 - 4*r01 + 3*r0);
     a = 1/(DT)^2*(2*r0 - 5*r01 + 4*r02 - 1*r03);
     da = 1/(DT)^3*(5/2*r0 - 9*r01 + 12*r02 - 7*r03 + 3/2*r04);
%      v = 1/2*(r02 - 4*r01 + 3*r0);
%      a = (2*r0 - 5*r01 + 4*r02 - 1*r03);
%      da = (5/2*r0 - 9*r01 + 12*r02 - 7*r03 + 3/2*r04);

     %t0 = v/norm(v);
     %t01 = -r02 + r0; %t01 = t01/norm(t01);
     %t02 = -r03 + r01; %t02 = t02/norm(t02);
     %t03 = -r04 - r02; %t03 = t03/norm(t03);
       
     %n0 = (t02 - 4*t01 + 3*t0);  
     %n0 = n0/norm(n0);
     %n01 = (-t02 + t0);  %n01 = n01/norm(n01);
     %n02 = (-t03 + t01);  %n02 = n02/norm(n02);
     
     %t = t0;
     %n = n0;
     
   else % Central difference
     v = (r11 - r01)/(2*DT);
     a = (r01 - 2*r0 + r11)/(DT)^2;
     da = 1/(DT)^3*(-1/2*r02 + 1*r01 - 1*r11 + 1/2*r12);
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

   end
     speed = norm(v);
     %if speed>127
     % warning('MyProgam:Limits','speed is large: %s',speed)
     %end
     accel = norm(a);
   
     dsdksi = sqrt(v'*v);
     d2sdksi2 = (v'*a)/dsdksi;
     d3sdksi3 = (((a'*a) + (v'*da))*dsdksi - (v'*a)*d2sdksi2)/(v'*v);
     
     drds = v/dsdksi;
     d2rds2 = 1/dsdksi^2*(a - d2sdksi2*drds);
     d3rds3 = 1/dsdksi^3*(da - d3sdksi3*drds - 3*dsdksi*d2sdksi2*d2rds2);
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
     rho = 1/kappa;
%      rho = speed^3/sqrt(speed^2*accel^2 - (drds'*d2rds2)^2);
     
     %drhodksi = -dsdksi*(d2rds2'*d3rds3)/(d2rds2'*d2rds2)^(3/2);
     %drhods = -(d2rds2'*d3rds3)/(d2rds2'*d2rds2)^(3/2);
%      rho3 = rho; rho4 = rho;
     tau = norm(1/kappa^2*(kappa*d3rds3 - dkappa*d2rds2) + drds*kappa);
     
     at = d2sdksi2;
     an = dsdksi^2/rho;
     
     T = drds;
     N = d2rds2/kappa;
     B = ((d3rds3*kappa - d2rds2*dkappa)/kappa^2 + drds*kappa)/tau; %B = B/norm(B);
     
     %rho = -N'*r0;
     rho = norm(v)^3/sqrt((norm(v)*norm(a))^2 - (v'*a)^2);
     differences=abs(rho-1/kappa);
     if differences>2.2088e-06
      warning('MyProgramm:lowPrecission','rho differs from 1/kappa by %s',differences)
     end
     
     %dN = (kappa*d3rds3 - dkappa*d2rds2)/kappa^2;

     %tau2 = B'*dN;
     
     ortCond2 = T'*N;
     ortCond3 = T'*B;
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
     rho2 = -N'*r0;

%    ortCond1 = abs((v'*v)./(r0'*a));
   ortCond1 = abs((r0'*a) + (v'*v));
   
   drddr = v'*a/norm(v)/norm(a);
   
   %cosmu = (r0'*a)/norm(a);
   cosmu = -((r0'*(r01 - 2*r0 + r11))/norm(r01 - 2*r0 + r11));
   cospsi = ((r11 - r01)'*(r11 - 2*r0 + r01))/(norm(r11 - r01)*norm(r11 - 2*r0 + r01));
   %sinpsi = sqrt(1 - (v'*a/norm(v)/norm(a)).^2);
   sinpsi = sqrt(1 - cospsi^2);
%    rho = abs(cosmu/sinpsi);
   R = r0;
   V = v;
   A = a;
   
end