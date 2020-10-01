function res = sortEigenValuesAndGetQuantities(model,sortType,plotfig,forcedeig)
   if nargin<1
        model = runEigenProblem();
        sortType = 'backwards';
%         sortType = 'forwards';
%         sortType = 'none';
        %plotfig = [1,2];
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
   %RHO3  = zeros(f,1);           %RHO3(k) = rhos3;
   %RHO4  = zeros(f,1);           %RHO4(k) = rhos4;
   TAU  = zeros(f,1);            %TAU(k) = taus;
   %TAU2  = zeros(f,1);           %TAU2(k) = tau2s;
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
   
   if kl==0
    warning('MyProgram:noStabilityLimit','assuming Stability-limit with last step')
    kl=numel(lambda0);
   end
   stability_limit = lambda0(kl) + lami;
   switch sortType
       case 'none'
           disp('No sorting of eigenvectors');
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
           r1 = reshape(r1,length(r1),1);
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

[v,a,speed,accel,at,an,rho,tau,~,rho2,drddr_,cosmu_,sinpsi_,ortCond1,ortCond2,ortCond3, ~, ~, B, R, V, A] ...
 = getQuantities(RS,dksi,epsilon); %DT=epsilon
           coplanar_ = B'*A/norm(A);
%            if rho>rhotry
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
      Orth1 = abs(OC1)<1e-4;
      if model.check==true
       S(isnan(S))=0;
       Orth2 = S<124;
       if ~all(Orth2)
        warning('MyProgam:Limits','speed is large: %f',max(S))
       end
       A0(isnan(A0))=0;
       Orth3 = A0<0.1;
       if ~all(Orth3)
        warning('MyProgam:Limits','total accerlation is large: %f',max(A0))
       end
       At(isnan(At))=0;
       Orth4 = abs(At)<0.035;
       if ~all(Orth4)
        warning('MyProgam:Limits','tangential accerlation is large: %f',max(At))
       end
       %       Orth5 = abs(An)<0.007;
       %       if ~all(Orth5)
       %        warning('MyProgam:Limits','normal accerlation is large: %f',max(An))
       %       end
       %Orth=logical(Orth1);
       Orth=logical((Orth1).*(Orth2).*(Orth3).*(Orth4));
      else
       Orth = true(length(Orth1),1);
      end
      
      res.lambdaorg = lambda';
      lambda = lambda(Orth);    res.lambda = lambda';
      V = V(:,Orth);            res.V = V;
      A = A(:,Orth);            res.A = A;
      S = S(Orth);              res.S = S;
      A0 = A0(Orth);            res.A0 = A0;
      At = At(Orth);            res.At = At;
      An = An(Orth);            res.An = An;
      res.RHO = RHO(Orth);
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
      res.cosmu = cosmu(Orth);
      Orthcosmu=[false;(res.cosmu(2:end)-res.cosmu(1:end-1))>0.6];
      res.cosmu(Orthcosmu)=NaN;
      sinpsi = sinpsi(Orth);    res.sinpsi = sinpsi;
      coplanar = coplanar(Orth); res.coplanar = coplanar;
      
      res.OC0 = OC0; res.OC1 = OC1; res.OC2 = OC2; res.OC3 = OC3;
      res.eigposition = eigposition;
      res.stability_limit = [kl,stability_limit];
      
      if isempty(plotfig)
       Disp('no plot')
      else
       MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
       colJK=MyColours{forcedeig};
       plotres(res,model,plotfig,colJK)
      end
end
