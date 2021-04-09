function rhorho3D = getQuantitiesStern(CosPhi,RHO2,lambda)

f=numel(lambda);
rhorho3D=NaN(f,1);
CosTheta=sqrt(1-RHO2.^2);
SinPhi=sqrt(1-CosPhi.^2);
for i = 3:f-1
   r01 = [CosTheta(i-1)*CosPhi(i-1);CosTheta(i-1)*SinPhi(i-1);RHO2(i-1)];
   r0  = [CosTheta( i )*CosPhi( i );CosTheta( i )*SinPhi( i );RHO2( i )];
   r11 = [CosTheta(i+1)*CosPhi(i+1);CosTheta(i+1)*SinPhi(i+1);RHO2(i+1)];
   DT1=lambda(i)-lambda(i-1);
   DT2=lambda(i+1)-lambda(i);
   DT=(DT1+DT2)/2;
   if abs(DT1-DT2)>eps(lambda(end))
    warning('MyProgam:Lambda','Lambda-steps differ')
    rhorho3D=NaN;
    return
   end
   


% Central difference
     v = (r11 - r01)/(2*DT);
     a = (r01 - 2*r0 + r11)/(DT)^2;

%      speed = norm(v);
   
     dsdksi = sqrt(v'*v);
%      diff=abs(speed-dsdksi);
%      limit=3e-14;
%      if ~isnan(diff) && diff>limit
%       if speed>124
%        warning('MyProgram:Boundary','Speed is large, returing partially NaN')
%        v=v*NaN;
%        a=a*NaN;
%        r0=NaN*r0;
%        dsdksi=NaN;
%       else
%        assert(diff<=limit,'norm(v) not equal sqrt(v*v) by %d',diff)
%       end
%      end
     d2sdksi2 = (v'*a)/dsdksi;
     drds = v/dsdksi;
     d2rds2 = 1/dsdksi^2*(a - d2sdksi2*drds);
     kappa = norm(d2rds2); %Gl19
     N = d2rds2/kappa;
     rhorho3D(i) = -N'*r0;

end
     
end