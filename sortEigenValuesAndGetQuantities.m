function res = sortEigenValuesAndGetQuantities(model,sortType,noplot,forcedeig)
   if nargin<1
        model = runEigenProblem();
        sortType = 'backwards';
%         sortType = 'forwards';
%         sortType = 'none';
        noplot = 0;
        forcedeig = [];
   end
   
   epsilon = model.lambda(2) - model.lambda(1);
   eigval = model.eigenvalues;
   eigvec = model.eigenvectors;
   lambda0 = model.lambda0;
   arclengths = model.arclengths;
   
   f = length(eigval);
   eigposition = zeros(length(lambda0),1);
   
   %% Computed Quantities:
   V    = zeros(size(eigvec{1},2),f); %V(:,k) = vs;
   A    = zeros(size(eigvec{1},2),f); %A(:,k) = as;     
   S    = zeros(f,1);            %S(k) = speeds;    
   A0   = zeros(f,1);            %A0(k) = a0s;    
   At   = zeros(f,1);            %At(k) = ats;      
   An   = zeros(f,1);            %An(k) = ans_;      
   RHO  = zeros(f,1);            %RHO(k) = rhos;
   RHO2  = zeros(f,1);           %RHO2(k) = rhos2;
   RHO3  = zeros(f,1);           %RHO3(k) = rhos3;
   RHO4  = zeros(f,1);           %RHO4(k) = rhos4;
   TAU  = zeros(f,1);            %TAU(k) = taus;
   TAU2  = zeros(f,1);           %TAU2(k) = tau2s;
   OC0  = zeros(f,1);            %OC0(k) = ortConds0;
   OC1  = zeros(f,1);            %OC1(k) = ortConds1;
   OC2  = zeros(f,1);            %OC2(k) = ortConds2;
   OC3  = zeros(f,1);            %OC3(k) = ortConds3;
   OCeig  = zeros(f,1);          %OCeig(k) = 1;
   LAM  = zeros(f,1);            %LAM(k) = lams;
   R    = zeros(size(eigvec{1},2),f); %R(:,k) = r0s;
   r1ri = zeros(f,1);
   rsri = zeros(f,1);
   drddr = zeros(f,1);
   cosmu = zeros(f,1);
   sinpsi = zeros(f,1);
   coplanar = zeros(f,1);
   
   % help quantities:
   POS = zeros(f,1);             %POS(k) = poslam;
   %%
   
   lambda = lambda0(1:f)';
   
   [kl,lami,~] = findStabilityLimit(eigval,lambda0);
   if lami<0
       kstability = kl - 1;
   else
       kstability = kl;
   end
   if kstability == 0
       kstability = 1;
   end
   
   stability_limit = lambda0(kl) + lami;
   switch sortType
       case 'none'
           disp('No sorting of eigenvectors');
       case 'forwards'
           k = 1; lami = 0; poslam = 1;
       case 'backwards'
           [k,lami,poslam] = findStabilityLimit(eigval,lambda0);
   end
   
   if ~strcmpi(sortType, 'none')
       r0try = eigvec{k}(5,:,poslam);  r0try = reshape(r0try,length(r0try),1);
       r1 = eigvec{1}(5,:,poslam);  r1 = reshape(r1,length(r1),1);
       rs = eigvec{kl}(5,:,poslam);  rs = reshape(rs,length(rs),1);
   else
       r0try = zeros(size(eigvec{1},2),1);
       if ~isempty(forcedeig)
           r1 = eigvec{1}(5,:,forcedeig);  r1 = reshape(r1,length(r1),1);
           rs = eigvec{kl}(5,:,forcedeig);  rs = reshape(rs,length(rs),1);
       else
           r1 = eigvec{1}(5,:,1);  r1 = reshape(r1,length(r1),1);
           rs = eigvec{kl}(5,:,1);  rs = reshape(rs,length(rs),1);
       end
   end
   
   for i = 1:f
       C0 = eigvec{i}(5,:,:); C0 = reshape(C0,size(C0,2),size(C0,3));
       
       if strcmpi(sortType, 'none')
           is0 = 1;
       else
           [val,is0] = compareEigenmodes(C0,r0try);
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
       
       rhotry = 0;
       for isi = 1:1%length(is)
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

           lami = eigval{i}(5,is0(isi));
           
           if ((lambda0(kstability)-2*epsilon)<=lambda0(i))&&(lambda0(i)<=lambda0(kstability))
               r11 = 0*r11; r12 = 0*r12; r13 = 0*r13; r14 = 0*r14;
           elseif (lambda0(kstability)<lambda0(i))&&(lambda0(i)<=(lambda0(kstability)+2*epsilon))
               r01 = 0*r01; r02 = 0*r02; r03 = 0*r03; r04 = 0*r04;
           elseif lambda0(i)<2*epsilon
               r01 = 0*r01; r02 = 0*r02; r03 = 0*r03; r04 = 0*r04;
           end
           RS = [r04, r03, r02, r01, r0, r11, r12, r13, r14];
           dksi = arclengths{i};
%            dksi = ones(size(dksi,1),size(dksi,2))*(lambda(2)-lambda(1));

           [v,a,speed,accel,at,an,rho,tau,ortCond0,rho2,drddr_,cosmu_,sinpsi_,ortCond1,ortCond2,ortCond3, T, N, B, R, V, A] = getQuantities(RS,dksi,epsilon);
           coplanar_ = B'*A/norm(A);
%            if rho>rhotry
               rhotry = rho;
               indi = 1;
       
               V(:,i) = v;
               A(:,i) = a;
               S(i) = speed;
               A0(i) = accel;
               At(i) = at;
               An(i) = an;
               RHO(i) = rho;
               RHO2(i) = rho2;
               TAU(i) = tau;
               OC0(i) = abs((v'*v)+(r0'*a));
               OC1(i) = ortCond1;
               OC2(i) = ortCond2;
               OC3(i) = ortCond3;
               OCeig(i) = abs(r0'*v)/norm(v);
               LAM(i) = lami;
        %        disp(indi);
               R(:,i) = indi*r0;
               POS(i) = is0(1);
               
               r1ri(i) = abs(r1'*r0);
               rsri(i) = abs(rs'*r0);
               drddr(i) = drddr_;
               cosmu(i) = cosmu_;
               sinpsi(i) = sinpsi_;
               coplanar(i) = coplanar_;
        %        r0try = indi*r0;
%            end
       end
   end
   
   %% Apply orthogonality conditions
      Orth = abs(OC1)<1e-4;
      Orth = true(length(Orth),1);
      
      res.lambdaorg = lambda';
      lambda = lambda(Orth);    res.lambda = lambda';
      V = V(:,Orth);            res.V = V;
      A = A(:,Orth);            res.A = A;
      S = S(Orth);              res.S = S;
      A0 = A0(Orth);            res.A0 = A0;
      At = At(Orth);            res.At = At;
      An = An(Orth);            res.An = An;
      RHO = RHO(Orth);          res.RHO = RHO;
      RHO2 = RHO2(Orth);        res.RHO2 = RHO2;
      TAU = TAU(Orth);          res.TAU = TAU;
%       OC1 = OC1(Orth);
      OCeig = OCeig(Orth);      res.OCeig = OCeig;
      LAM = LAM(Orth);          res.LAM = LAM;
      R = R(:,Orth);            res.R = R;
      POS = POS(Orth);          res.POS = POS;
      r1ri = r1ri(Orth);        res.r1ri = r1ri;
      rsri = rsri(Orth);        res.rsri = rsri;
      drddr = drddr(Orth);      res.drddr = drddr;
      cosmu = cosmu(Orth);      res.cosmu = cosmu;
      sinpsi = sinpsi(Orth);    res.sinpsi = sinpsi;
      coplanar = coplanar(Orth); res.coplanar = coplanar;
      
      res.OC0 = OC0; res.OC1 = OC1; res.OC2 = OC2; res.OC3 = OC3;
      res.eigposition = eigposition;
      res.stability_limit = [kl,stability_limit];
      
      if noplot == 0
          plotres(res)
      end
end
function plotres(res)
      lambda1 = res.lambdaorg;
      lambda = res.lambda;
      V = res.V;
      A = res.A;
      S = res.S;
      A0 = res.A0;
      At = res.At;
      An = res.An;
      RHO = res.RHO;
      RHO2 = res.RHO2;
      TAU = res.TAU;
      cosmu = res.cosmu;
      sinpsi = res.sinpsi;
%       OC1 = OC1(Orth);
      OCeig = res.OCeig;
      LAM = res.LAM;
      R = res.R;
      POS = res.POS;
      
      OC0 = res.OC0;
      OC1 = res.OC1; 
      OC2 = res.OC2; 
      OC3 = res.OC3;
      
      coplanar = res.coplanar;

      colo = [rand(), rand(), rand()];  
      colo = [0 0 0];
      
      kl = res.stability_limit(1);
      stability_limit = res.stability_limit(2);
      
   if nargin<2
       maxlam = max(lambda);
       maxlam = 1.0;
       lambda = lambda/maxlam;
       figure(2); 
       hold on
       plot(lambda(2:end),RHO(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',colo);
       
%        plot(lambda(2:end),RHO2(2:end),'LineStyle','none','Marker','.','Color',colo);
%        plot(lambda,RHO3,'LineStyle','none','Marker','o','Color',colo);
%        plot(lambda,RHO4,'LineStyle','--','Marker','x','Color',colo);
       xlabel('lambda');
       ylabel('rho');
       title('rho');
       bbb = gca();
       bbb.YLim = [0.0,1.1];
       col = bbb.Children(1).Color;
       
%        plot(lambda(kl),RHO(kl),'LineStyle','none','Marker','x','LineWidth',1.5,'Color',col);
       grid on

       figure(1);
       hold on
       lambda = reshape(lambda,length(lambda),1); LAM = reshape(LAM,length(lambda),1);
       LAM = LAM/maxlam;
%        plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none');%,'Color',colo);
%        plot(lambda(2:end), LAM(2:end),'LineStyle','--','Marker','none');%,'Color',colo);
       aaa = gca();
       if size(aaa.Children,1)<2
          plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
          plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
       end
       plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
       xlabel('Lambda');
%        ylabel('Lambda + chi | or | chi');
       ylabel('chi');
       grid on
%        legend('chi','eigenvalue','lambda','zero')

       figure(3)
       hold on
       plot(lambda(2:end),S(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
       xlabel('lambda');
       ylabel('velocity');
       title('velocity');
       grid on
       nopeaks = (S(S<10*mean(S(:))));
       bbb = gca();
       bbb.YLim = [0.0,max(nopeaks)];

       figure(4)
       hold on
       At = real(At);
       plot(lambda(2:end),At(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
       xlabel('lambda');
       ylabel('tangential acceleration');
       title('tangential acceleration');
       grid on
       nopeaks = (At(abs(At)<10*mean(abs(At(:)))));
       bbb = gca();
       bbb.YLim = [min(nopeaks),max(nopeaks)];

       figure(5)
       hold on
       plot(lambda(2:end),An(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
       xlabel('lambda');
       ylabel('normal acceleration');
       title('normal acceleration');
       grid on
       nopeaks = (An(An<10*mean(An(:))));
       bbb = gca();
       bbb.YLim = [0.0,max(nopeaks)];
       
       figure(6)
       hold on
       plot(lambda(2:end),A0(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
       xlabel('lambda');
       ylabel('total acceleration');
       title('total acceleration');
       grid on
       nopeaks = (A0(A0<10*mean(A0(:))));
       bbb = gca();
       bbb.YLim = [0.0,max(nopeaks)];

       figure(7)
       hold on
       plot(lambda(2:end),abs(TAU(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);%,'Color',colo);
%       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o','Color',colo);
       xlabel('lambda');
       ylabel('torque');
       title('torque');
       grid on
       nopeaks = (TAU(TAU<10*mean(TAU(:))));
       bbb = gca();
       bbb.YLim = [0.0,max(nopeaks)];
       
       figure(8)
       hold on
       plot(lambda(2:end),abs(sinpsi(2:end)),'LineStyle','-','Marker','x','LineWidth',1.5,'Color',col);
       plot(lambda(2:end),abs(cosmu(2:end)),'LineStyle','-','Marker','o','LineWidth',1.5,'Color',col);
       xlabel('lambda');
       ylabel('absolute values of: sin(psi) and cos(mu)');
       title('absolute values of: sin(psi) and cos(mu)');
       grid on
       
% % %        figure(9)
% % %        hold on
% % %        plot(lambda(2:end),asin(sinpsi(2:end))*180/pi,'LineStyle','-','LineWidth',1.5,'Marker','x','Color',col);
% % %        plot(lambda(2:end),acos(cosmu(2:end))*180/pi,'LineStyle','--','LineWidth',1.5,'Marker','o','Color',col);
% % % %        plot(lambda(2:end),(asin(sinpsi(2:end))+acos(cosmu(2:end)))*180/pi,'LineStyle','-','Marker','.','Color',[1,0,0]);
% % % %        plot(lambda(2:end),abs(asin(sinpsi(2:end))-acos(cosmu(2:end)))*180/pi,'LineStyle','--','Marker','.','Color',[1,0,0]);
% % %        xlabel('Lambda');
% % %        ylabel('angle');
% % % %        legend('psi','nu','psi+nu','abs(psi-nu)')
% % %        grid on
% % % %        nopeaks = (TAU(TAU<10*mean(TAU(:))));
% % %        bbb = gca();
% % %        maxy = max(bbb.YTick);
% % %        bbb.YTick = 0.0:45:maxy;
       
% % %        figure(10)
% % %        if (asin(sinpsi(2))*180/pi>=0)&&(acos(cosmu(2))*180/pi<=90)
% % %            psi1 = asin(sinpsi)*180/pi;
% % %            nu1 = acos(cosmu)*180/pi;
% % %            psi2 = asin(sinpsi)*180/pi;
% % %            nu2 = 360 - acos(cosmu)*180/pi;
% % %            psi3 = 180 - asin(sinpsi)*180/pi;
% % %            nu3 = acos(cosmu)*180/pi;
% % %            psi4 = 180 - asin(sinpsi)*180/pi;
% % %            nu4 = 360 - acos(cosmu)*180/pi;
% % %        elseif (asin(sinpsi(2))*180/pi>=0)&&(acos(cosmu(2))*180/pi>=90)
% % %            psi1 = asin(sinpsi)*180/pi;
% % %            nu1 = acos(cosmu)*180/pi;
% % %            psi2 = asin(sinpsi)*180/pi;
% % %            nu2 = 360 - acos(cosmu)*180/pi;
% % %            psi3 = 180 - asin(sinpsi)*180/pi;
% % %            nu3 = acos(cosmu)*180/pi;
% % %            psi4 = 180 - asin(sinpsi)*180/pi;
% % %            nu4 = 360 - acos(cosmu)*180/pi;
% % %        elseif (asin(sinpsi(2))*180/pi<=0)&&(acos(cosmu(2))*180/pi<=90)
% % %            psi1 = 180 - asin(sinpsi)*180/pi;
% % %            nu1 = acos(cosmu)*180/pi;
% % %            psi2 = 180 - asin(sinpsi)*180/pi;
% % %            nu2 = 360 - acos(cosmu)*180/pi;
% % %            psi3 = 360 + asin(sinpsi)*180/pi;
% % %            nu3 = acos(cosmu)*180/pi;
% % %            psi4 = 360 + asin(sinpsi)*180/pi;
% % %            nu4 = 360 - acos(cosmu)*180/pi;
% % %        elseif (asin(sinpsi(2))*180/pi<=0)&&(acos(cosmu(2))*180/pi>=90)
% % %            psi1 = 180 - asin(sinpsi)*180/pi;
% % %            nu1 = acos(cosmu)*180/pi;
% % %            psi2 = 180 - asin(sinpsi)*180/pi;
% % %            nu2 = 360 - acos(cosmu)*180/pi;
% % %            psi3 = 360 + asin(sinpsi)*180/pi;
% % %            nu3 = acos(cosmu)*180/pi;
% % %            psi4 = 360 + asin(sinpsi)*180/pi;
% % %            nu4 = 360 - acos(cosmu)*180/pi;
% % %        end
% % %        subplot(2,2,1)
% % %        hold on
% % %        plot(lambda(2:end),psi1(2:end),'LineStyle','-','Marker','x','Color',colo);
% % %        plot(lambda(2:end),nu1(2:end),'LineStyle','--','Marker','o','Color',colo);
% % %        plot(lambda(2:end),psi1(2:end)+nu1(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
% % %        bbb = gca();
% % %        bbb.YLim = [0,360];
% % %        bbb.YTick = 0.0:90:360;
% % %        grid on
% % %        
% % %        subplot(2,2,2)
% % %        hold on
% % %        plot(lambda(2:end),psi2(2:end),'LineStyle','-','Marker','x','Color',colo);
% % %        plot(lambda(2:end),nu2(2:end),'LineStyle','--','Marker','o','Color',colo);
% % %        plot(lambda(2:end),psi2(2:end)+nu2(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
% % %        bbb = gca();
% % %        bbb.YLim = [0,360];
% % %        bbb.YTick = 0.0:90:360;
% % %        grid on
% % %        
% % %        subplot(2,2,3)
% % %        hold on
% % %        plot(lambda(2:end),psi3(2:end),'LineStyle','-','Marker','x','Color',colo);
% % %        plot(lambda(2:end),nu3(2:end),'LineStyle','--','Marker','o','Color',colo);
% % %        plot(lambda(2:end),psi3(2:end)+nu3(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
% % %        bbb = gca();
% % %        bbb.YLim = [0,360];
% % %        bbb.YTick = 0.0:90:360;
% % %        grid on
% % %        
% % %        subplot(2,2,4)
% % %        hold on
% % %        plot(lambda(2:end),psi4(2:end),'LineStyle','-','Marker','x','Color',colo);
% % %        plot(lambda(2:end),nu4(2:end),'LineStyle','--','Marker','o','Color',colo);
% % %        plot(lambda(2:end),psi4(2:end)+nu4(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
% % %        bbb = gca();
% % %        bbb.YLim = [0,360];
% % %        bbb.YTick = 0.0:90:360;
% % %        grid on
% % %        
% % %        figure(3000)
% % %        hold on
% % %        plot(lambda1,abs(OC0),'LineStyle','--','Marker','s','Color',colo);
% % %        hold on
% % %        plot(lambda1,abs(OC1),'LineStyle','--','Marker','o','Color',colo);
% % %        plot(lambda1,abs(OC2),'LineStyle','--','Marker','.','Color',colo);
% % %        plot(lambda1,abs(OC3),'LineStyle','--','Marker','^','Color',colo);
% % %        plot(lambda,1-OCeig,'LineStyle','-','Marker','x','Color',colo);
% % %        plotaxes = gca();
% % %        plotaxes.YScale = 'log';
% % %        xlabel('Lambda');
% % %        ylabel('Orthogonality checks');
% % %     %    legend({'n*t','1 - (r_i * r_{i+1})'})
% % %        legend({'n*t','T*N','T*B','N*B','1 - (r_i * r_s)'})
% % %        grid on
% % %        hold off
% % %        
% % %        figure(3001)
% % %        hold on
% % %        plot(lambda(2:end),OCeig(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);
% % %        plot(lambda(2:end),OC0(2:end),'LineStyle','--','Marker','none','LineWidth',1.5,'Color',col);
% % %        plotaxes = gca();
% % % %        plotaxes.YScale = 'log';
% % %        xlabel('Lambda');
% % %     %    legend({'n*t','1 - (r_i * r_{i+1})'})
% % % %        legend({'r * r'''})
% % %        grid on

   end
end