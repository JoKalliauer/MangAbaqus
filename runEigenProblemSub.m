function model = runEigenProblemSub(modelprops,model,Displ,Kts,Kg,matches,wbrEP)
 
 % [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
 
 %Displ = NodalResults2Displ(Nres);
 if numel(Displ)<1
  displacementsenable=false;
 else
  displacementsenable=true;
 end
 %Kg = EigRes;
 lambda0=model.lambda0';
 eigval = cell(length(lambda0),1);
 eigvec = cell(length(lambda0),1);
 displacements = cell(length(lambda0),1);
 darclengths = cell(length(lambda0),1);
 arclengthJK = cell(length(lambda0),1);
 arclengthHM = cell(length(lambda0),1);
 StiffMtxs = cell(length(lambda0),3);
 DetKtx = NaN(length(lambda0),1);
 %load0 = NaN(length(lambda0),1);
 eigvecDRH = cell(length(lambda0),1);% DRH...Displacement,Rotation,Hybrid(splitted)
 
%  matches = NaN(0);
%  n = 0;
 %m = 0;
%  iend=min((size(num,1)-1),size(fulllambda,1));
%  for i = 1:iend
%   if ~isempty(find(lambda0==fulllambda(i), 1))
%    n = n+1;
%    matches(n) = i;
%   end
%  end
%  if matches(1)~=1
%   matches = [1,matches] %#ok<NOPRT>
%  end
 
 %Energy = zeros(length(lambda0),3);
 %Energy(:,1) = lambda0;
 mintest=max(min(size(Kts,1),size(Displ,1)),3);
%  if max(matches)+2>mintest
%   %warning('MyProgram:Abaqus','Abaqus might exited with error (step %d, l=%f) befor last lamdba (l=%f)',size(Kts,1),model.fulllambda(min(size(Kts,1),numel(model.fulllambda))),max(model.fulllambda))
%   matches(matches+2>mintest)=[];
%  end
 
 %% solve EigvalueProblem
 Kt0_0 = Kts{matches(1),2};
 oldsizeKt0=size(Kts{matches(1),2},1);
 ru36 = diag(Kt0_0==1e36);% remove boundary conditions
 ruinf = diag(Kt0_0==inf);
 ru=logical(ru36+ruinf);
 BC=find(ru);
 assert(all(BC==model.BC(:,1)),'ru does not agree with BC')
 Kt0_0(ru,:) = []; Kt0_0(:,ru) = [];
 newsizeKt0=size(Kt0_0,1);
 %[~ ,~ ,numofeigs0 ] = solveCLEforMinEigNew(Kt0_0,Kt0_0,Kg,Kt0_0,modelprops.typeofanalysis,matches(1),NaN,modelprops);
 Kt11= Kts{matches(1)+1,2};
 Kt11(ru,:) = []; Kt11(:,ru) = [];
 Ktprim0_0 = 1/(modelprops.epsilon)*(Kt11 - Kt0_0);
 [r0t0 ,~ ,numofeigs0 ] = solveCLEforMinEigNew(Kt11,Ktprim0_0,Kg,Kt0_0,modelprops.typeofanalysis,matches(1)+1,model,modelprops,Ktprim0_0);
 numofeigs=min([modelprops.numofeigs,newsizeKt0,numofeigs0]);
 newra = transpose(1:newsizeKt0);
 fullEV=NaN(numofeigs,size(model.fulllambda,1));
 RR0 = NaN(newsizeKt0,numofeigs);
 if numel(Displ)>0
  sizeDisp=size(Displ{1},1);
  if size(Displ{1},1)~=oldsizeKt0
   warning('MyProgram:Disp','size(Displ{1},1) has a different size than oldsizeKt0')
  end
 else
  warning('MyProgram:Disp','Displ from *.dat are empty, maybe Abaqus-Error')
  sizeDisp=min(oldsizeKt0);
 end
 displacements_=NaN(sizeDisp,9);
 
 model.N0=size(Kt0_0,1);
 
 if strcmp(modelprops.elementtype(1:2),'B2')
  model.RestrictedDOFs=3:5;
  model.Dim=2;
 else
  model.RestrictedDOFs=[];
  model.Dim=3;
 end
 
 R = NaN(7,size(r0t0,1),size(r0t0,2)); % dl x DoF x NrEigs
 R_DRsize=[5,model.dofpNode,size(model.Nodes,1),size(r0t0,2)]; % dl x DoFpNode x Nodes x NrEigs
 R_DRHsize=R_DRsize+[0 -numel(model.RestrictedDOFs) model.inDOF(2)-model.inDOF(1)+1 0]; % dl x DoFpNode(reducedRestricted) x Nodes(inkl.Hyb) x NrEigs
 R_DRH = NaN(R_DRHsize); % dl x DoFpNode x Nodes x NrEigs
 R_DRHsizeIn=R_DRHsize(2:4);%+[0 inNr 0]; %DoFpNode x Nodes(inkl.Hyb) x NrEigs
 ReducedHybridDofa=model.inDOF(3)+1:R_DRHsizeIn(1);
 ReducedHybridDofb=model.inDOF(4)+1:R_DRHsizeIn(1);
 nA=false(R_DRHsizeIn(1:2));
 nA(ReducedHybridDofa,model.inDOF(1):2:model.inDOF(2))=true;
 if ~isnan(ReducedHybridDofb)
  nA(ReducedHybridDofb,model.inDOF(1)+1:2:model.inDOF(2))=true;
 end
 if sum(strcmp(fieldnames(model), 'JC')) == 0
  model.JC=[];
 else
  for i1=1:size(model.JC,1)
   nA(model.JC(i1,2),model.JC(i1,1))=true;%only works if model.JC(i,2) is smaller than model.RestrictedDOFs
  end
 end
 if sum(strcmp(fieldnames(modelprops), 'sigma')) == 0
  modelprops.sigma=0;
 end
 
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblem');end

 for i = 1:length(matches)
  if usejava('jvm'); waitbar(i/length(matches),wbrEP,'runEigenProblem EigvalueProblem');end
  %disp('Lambda:');
  %disp(lambda(matches(i)));
  %Lambda=fulllambda(matches(i)) %#ok<NASGU,NOPRT>
  %disp('-------------');
  %     Energy(i,2) = membrane(matches(i));
  %     Energy(i,3) = nonmembrane(matches(i));
 %for i
  if i==0

   
  else % if i~= 1
   % -4.
   Zeile=matches(i)-4;
   if Zeile>0
    if Zeile>1
     dksi04 = sqrt((Displ{Zeile}-Displ{matches(i)-5})'*(Displ{Zeile}-Displ{matches(i)-5}));
     displacements_(:,1) = Displ{matches(i)-5};
    else
     dksi04 = NaN;
     displacements_(:,1) = NaN*Displ{1};
    end
    %Kt04 = Kts{Zeile,2};
    dksi03 = sqrt((Displ{matches(i)-3}-Displ{matches(i)-4})'*(Displ{matches(i)-3}-Displ{matches(i)-4}));
    displacements_(:,2) = Displ{matches(i)-4};
   else
    %Kt04 = 0*Kts{1,2};
    dksi04 = NaN;
    dksi03 = NaN;
    if displacementsenable
     displacements_(:,1) = NaN*Displ{1};
     displacements_(:,2) = NaN*Displ{1};
    end
   end
   %Kt04(ru,:) = []; Kt04(:,ru) = [];
   % -3.
   if Zeile>-1
    Kt03 = Kts{matches(i)-3,2};
    dksi02 = sqrt((Displ{matches(i)-2}-Displ{matches(i)-3})'*(Displ{matches(i)-2}-Displ{matches(i)-3}));
    displacements_(:,3) = Displ{matches(i)-3};
   else
    Kt03 = 0*Kts{1,2};
    dksi02 = NaN;
    displacements_(:,3) = NaN;
   end
   Kt03(ru,:) = []; Kt03(:,ru) = [];
   % -2.
   if Zeile>-2
    Kt02 = Kts{matches(i)-2,2};
    dksi01 = sqrt((Displ{matches(i)-1}-Displ{matches(i)-2})'*(Displ{matches(i)-1}-Displ{matches(i)-2}));
    displacements_(:,4) = Displ{matches(i)-2};
   else
    Kt02 = 0*Kts{1,2};
    dksi01 = NaN;
    displacements_(:,4) = NaN;
   end
   Kt02(ru,:) = []; Kt02(:,ru) = [];
   % -1.
   if Zeile>-3
    Kt01 = Kts{matches(i)-1,2};
    dksi11 = sqrt((Displ{matches(i)}-Displ{matches(i)-1})'*(Displ{matches(i)}-Displ{matches(i)-1}));
    displacements_(:,5) = Displ{matches(i)-1};
   else
    Kt01 = 0*Kts{1,2};
    dksi11 = NaN;
    displacements_(:,5) = NaN;
   end
   Kt01(ru,:) = []; Kt01(:,ru) = [];
   % 0.
   KT = Kts{matches(i),2}; KT(ru,:) = []; KT(:,ru) = [];
   % 1.
   Kt11 = Kts{matches(i)+1,2};
   if numel(Kt11)>0
    Kt11(ru,:) = []; Kt11(:,ru) = [];
   else
    Kt11=0*KT;
   end
   % 2.
   if displacementsenable
    dksi12 = sqrt((Displ{matches(i)+1}-Displ{matches(i)})'*(Displ{matches(i)+1}-Displ{matches(i)}));
   else
    dksi12 = NaN;
   end

   % 3.
   %dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));

   % 4.
   if displacementsenable
   displacements_(:,6) = Displ{matches(i)-0};
   displacements_(:,7) = Displ{matches(i)+1};
   if size(Displ,1)<matches(i)+2
    displacements_(:,8)=NaN*displacements_(:,7);
   else
    displacements_(:,8) = Displ{matches(i)+2};
   end
   %displacements_(:,9) = Displ{matches(i)+3};
   end
   
   %Ktprim03 = 1/(2*epsil)*(Kt02 - Kt04);
   Ktprim02 = 1/(2*modelprops.epsilon)*(Kt01 - Kt03);
   Ktprim01 = 1/(2*modelprops.epsilon)*(KT - Kt02);
   Ktprim0 = 1/(2*modelprops.epsilon)*(Kt11 - Kt01);



   if size(Kts,1)<matches(i)+2
    Ktprim11=0*Ktprim0;
   else
    Kt12 = Kts{matches(i)+2,2};
    if numel(Kt12)>0
     Kt12(ru,:) = []; Kt12(:,ru) = [];
    else
     Kt12=0*Kt11;
    end
    Ktprim11 = 1/(2*modelprops.epsilon)*(Kt12 - KT);
   end
   
   if size(Kts,1)<matches(i)+3
    Ktprim12=0*Ktprim11;
   else
    Kt13 = Kts{matches(i)+3,2};
    if numel(Kt13)>0
     Kt13(ru,:) = [];
     Kt13(:,ru) = [];
     Ktprim12 = 1/(2*modelprops.epsilon)*(Kt13 - Kt11);
    else
     Ktprim12=0*Ktprim11;
    end
   end
   
   if strcmp(modelprops.typeofanalysis,'KNL2') || mintest<matches(i)+4
    %Kt14=0*Kt13;
    %Ktprim13=0*Ktprim12;
    %dksi14=NaN;
    if displacementsenable; displacements_(:,9)=NaN*displacements_(:,8);end
   else
    %Kt14 = Kts{matches(i)+4,2};
    %Kt14(ru,:) = [];
    %Kt14(:,ru) = [];
    %Ktprim13 = 1/(2*epsil)*(Kt14 - Kt12);
    %Ktprim14 = 1/epsil*(3/2*Kt14 - 2*Kt13 + 1/2*Kt12);
    %dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
    displacements_(:,9) = Displ{matches(i)+3};
   end
   
   dksi = [dksi04, dksi03, dksi02, dksi01, dksi11,dksi12];
   darclengths{i} = dksi;
   
   displacements_(abs(displacements_)<=1e-31)=0;%remove numeric issues close to zero
   displacements{i} = displacements_;
   arclengthJK{i}=mean(sqrt(displacements_(:,5).*displacements_(:,5)));
   arclengthHM{i}=sqrt(sum(displacements_(:,5).*displacements_(:,5)));
   
   if i>1
    if matches(i)==matches(i-1)+1
     foloworder=true;
    else
     foloworder=false;
    end
   else
    foloworder=false;
   end
   
   %dKtprim0 = 1/epsil*(Kt01 -2*Kt0 +Kt11);
   if strcmp(modelprops.typeofanalysis,'Kg')
    [r0t ,eigval0_ ] = solveCLEforMinEigNew(KT ,Ktprim0 ,Kg,Kt0_0,modelprops.typeofanalysis,matches(i),model,modelprops);
    R = NaN(7,size(r0t,1),size(r0t,2));
    R(5,:,:) = r0t;
    EV = NaN(7,length(eigval0_));
    EV(5,:) = eigval0_;
   else
    if foloworder==true
     %r03=r02_
     r02_=r01_;
     r01_=r0_;
     r0_=r11_;
     r11_=r12_;
     eigval02_=eigval01_;
     eigval01_=eigval0_;
     eigval0_=eigval11_;
     eigval11_=eigval12_;
    else
     %[r03_,eigval03_] = solveCLEforMinEigNew(Kt03,Ktprim03,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-3,NaN,modelprops);
     [r02_,eigval02_] = solveCLEforMinEigNew(Kt02,Ktprim02,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-2,model,modelprops,Ktprim0_0);
     [r01_,eigval01_] = solveCLEforMinEigNew(Kt01,Ktprim01,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)-1,model,modelprops,Ktprim0_0);
     [r0_ ,eigval0_ ] = solveCLEforMinEigNew(KT ,Ktprim0 ,Kg,Kt0_0,modelprops.typeofanalysis,matches(i),model,modelprops,Ktprim0_0); % Kt=Kt0; iter=matches(i); typeofanal=modelprops.typeofanalysis
     [r11_,eigval11_] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+1,model,modelprops,Ktprim0_0);
    end
    [r12_,eigval12_] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+2,model,modelprops,Ktprim0_0);

%    if numel(Kt13)>0
%     [r13_,eigval13_] = solveCLEforMinEigNew(Kt13,Ktprim13,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+3,NaN,modelprops);
%     [r14_,eigval14_,numofeigs14] = solveCLEforMinEigNew(Kt14,Ktprim14,Kg,Kt0_0,modelprops.typeofanalysis,matches(i)+4,NaN,modelprops);
%    else
%     r13_=NaN*r12_;
%     eigval13_=NaN*eigval12_;
%    end
   
   
   
   if size(RR0,1)~=length(r0_)
    %r03t = RR0;
    r02t = RR0; r01t = RR0; r0t = RR0; r11t = RR0; r12t = RR0; %r13t = RR0; % r14t = RR0;
    %r03t(newra,:) = r03_;
    r02t(newra,:) = r02_;
    r01t(newra,:) = r01_;
    r0t(newra,:) = r0_;
    r11t(newra,:) = r11_;
    r12t(newra,:) = r12_;
    %r13t(newra,:) = r13_;
    %r14t(newra,:) = r14_;
   else
    %r03t = r03_;%/norm(r03_);
    r02t = r02_;%/norm(r02_);
    r01t = r01_;%/norm(r01_);
    r0t = r0_;%/norm(r0_);
    r11t = r11_;%/norm(r11_);
    r12t = r12_;%/norm(r12_);
    %r13t = r13_;%/norm(r13_);
    %r14t = r14_/norm(r14_);
   end
   
   
   for k = 1:size(r0_,2)
    if r01t(:,k)'*r0t(:,k)<0;  r01t(:,k)  = -r01t(:,k);  end
    if r02t(:,k)'*r01t(:,k)<0; r02t(:,k) = -r02t(:,k); end
    %if r03t(:,k)'*r02t(:,k)<0; r03t(:,k) = -r03t(:,k); end
    if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
    if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
    %if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end
    %if r13t(:,k)'*r14t(:,k)<0; r14t(:,k) = -r14t(:,k); end
   end
   
   %R(2,:,:) = r03t;
   R(3,:,:) = r02t;
   R(4,:,:) = r01t;
   R(5,:,:) = r0t;
   R(6,:,:) = r11t;
   R(7,:,:) = r12t;
   %R(8,:,:) = r13t;
   %R(9,:,:) = r14t;
   R_DRH(1,:,:,:)=RealEV(r02t,BC,R_DRHsize,nA);
   R_DRH(2,:,:,:)=RealEV(r01t,BC,R_DRHsize,nA);
   R_DRH(3,:,:,:)=RealEV(r0t ,BC,R_DRHsize,nA);
   R_DRH(4,:,:,:)=RealEV(r11t,BC,R_DRHsize,nA);
   R_DRH(5,:,:,:)=RealEV(r12t,BC,R_DRHsize,nA);
   if max(abs(imag(R)))>eps(0)
    warning('MyProgram:Complex','R is komplex');
   end
   
   EV = NaN(7,length(eigval0_)); %overwrite everything with NaNs
   %EV(2,:) = eigval03_;
   EV(3,:) = eigval02_;
   EV(4,:) = eigval01_;
   EV(5,:) = eigval0_;
   EV(6,:) = eigval11_;
   EV(7,:) = eigval12_;
   %EV(8,:) = eigval13_;
   %EV(9,:) = eigval14_;
   end
   if strcmp(modelprops.typeofanalysis,'KNL2')
    modelprops.sigma=min(modelprops.sigma,min(real(EV(:))));
   end
   if modelprops.followsigma==true && modelprops.numofeigs==1
    modelprops.sigma=eigval0_;
   end

   
  
  end
  
  if exist('eigval11_','var')
   if matches(i)<4
    %first=max(matches(i)-3,1);
    fullEV(:,matches(i):matches(i)+2)=[eigval0_,eigval11_,eigval12_];
    %fullLAMDA(matches(i):matches(i)+3)=fulllambda(matches(i):matches(i)+3);
    if matches(i)>=2
     fullEV(:,matches(i)-1)=eigval01_;
     %fullLAMDA(matches(i):matches(i)+3)
     if matches(i)>=3
      fullEV(:,matches(i)-2)=eigval02_;
     end
    end
   else
    fullEV(:,matches(i)-2:matches(i)+2)=[eigval02_,eigval01_,eigval0_,eigval11_,eigval12_];
    %fullLAMDA(matches(i)-3:matches(i)+3)=matches(i)-3:matches(i)+3;
   end % matches(i)<4
  else
   fullEV(:,matches(i))=eigval0_(1);
  end
  

  eigval{i} = EV;
  eigvec{i} = R;
  eigvecDRH{i}=R_DRH;% DRH...[Displacement,Rotation,Hybrid](splitted)
  StiffMtxs{i,1} = KT;
  StiffMtxs{i,2} = Ktprim0;
  if i==1
   %    KTmult0=5*10^-11;
   KTmult0=1.1e-11;
   KTrel=KT*KTmult0;
   res=det(KTrel);
   relative=false;
   if (res==0 && ~strcmp(modelprops.elementtype(end),'H')) || isinf(res)
    Kt0=KT;
    Kt0inv=inv(Kt0);
    relative=true;
%     n=87;
%     passt=false;
%     while passt==false
%      res=1;
%      %n=1;
%      while res>0 && ~isinf(res)
%       n=n+1;
%       res=det(Kt0(1:n,1:n)*KTmult0);
%      end
%      m=n-1
%      res=det(Kt0(1:m,1:m)*KTmult0)
%      if isinf(res)
%       warning('MyPrgm:Nan','res is inf')
%       %res=det(Kt0(1:m,1:m)*KTmultAlt);
%      end
%      if res<1e10
%       res=res*1e50 %1e100->2813
% %1e50: 1931->2189->2486
%      end
%      Nenner=nthroot(res,m);
%      %KTmultAlt=KTmult0;
%      KTmultNeu=(KTmult0/Nenner+KTmult0)/2
%      KTmult=KTmultNeu;
%      KTmult0=KTmult;
%      detKtmult=det(Kt0*KTmult);
%      if detKtmult~=0 && ~isinf(detKtmult)
%       passt=true;
%      end
%     end
   else
    KTmult=KTmult0;
   end
 end %i=1
 if relative==true
  if numel(KT)<=1E5
   DetKtx(i)=det(Kt0inv*KT); %#ok<MINV>
  else
   DetKtx(i)=NaN;
  end
 else
  KTrel=KT*KTmult;
  if isinf(res) || res==0
   DetKtx(i)=res;
  else
   DetKtx(i)=det(KTrel);
  end
 end
  %load(i)=model.fullload(matches(i));
  
  %StiffMtxs{i,3} = dKtprim0;
  
  %    [eq1,eq2] = solcheck(Kt0, typeofanal, Ktprim0, dKtprim0, lambda(matches(i)), epsil, EV, R, i, model);
  %     for p = 1:10
  %         r = R(5,:,p); r = r(:);
  %         r(ru) = [];
  %         L = EV(5,p);
  %         disp(norm((Kt0 + L*Ktprim0)*r));
  %     end
 end%for i = 1:length(matches)
 if usejava('jvm'); waitbar(1,wbrEP,'runEigenProblem EigvalueProblem finsih');end
 
 model.eigenvalues = eigval;
 model.eigenvectors = eigvec;
 model.eigvecDRH=eigvecDRH;% DRH...Displacement,Rotation,Hybrid(splitted)
 %model.energy = Energy;
 model.arclengths = darclengths;
 model.arclengthurJK = arclengthJK;
 model.arclengthurHM = arclengthHM;
 model.displacements = displacements;
 model.lambda0 = lambda0';
 model.load0=[0;model.load];
 model.stiffnessMatrices = StiffMtxs;
 model.DetKtx=DetKtx;
 %model.N0=size(Kt0_0,1);
 model.N=size(KT,2);
 assert(model.N0==model.N,'size of Kt0 changes or is not symmetric')
 %model.fullLAMDA=fullLAMDA;
 model.fullEV=fullEV;
 model.numofeigs=numofeigs;
 

 
end %fucntion

