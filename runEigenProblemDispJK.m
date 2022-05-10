function model = runEigenProblemDispJK(modelprops,model,Displ,~,~,matches,wbrEP)
 
 % [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
 
 %Displ = NodalResults2Displ(Nres);
 if numel(Displ)<1
  displacementsenable=false;
 elseif  numel(Displ)==2
  warning('MyPrgm:Disp:DiplOnly2','Displ only has 2 entries')
  displacementsenable=false;
 else
  displacementsenable=true;
 end
 %Kg = EigRes;
 lambda0=model.lambda0';

 displacements = cell(length(matches),1);
 darclengths = cell(length(lambda0),1);
 arclengthJK = cell(length(lambda0),1);
 lengthMaxJK = cell(length(lambda0),1);
 wMaxJK = cell(length(lambda0),1);
 vMaxJK = cell(length(lambda0),1);
 uMaxJK = cell(length(lambda0),1);
 arclengthHM = cell(length(lambda0),1);
 dxidl = cell(length(lambda0),1);
 dUdl = cell(length(lambda0),1);
 dwdl = cell(length(lambda0),1);
 dvdl = cell(length(lambda0),1);
 dudl = cell(length(lambda0),1);
 sDiffJK = NaN(length(lambda0),1);
 idx1= NaN(length(lambda0),1);
 idx2= NaN(length(lambda0),1);
 duEdl= NaN(length(lambda0),1);
 
 NrNo=size(model.Nodes,1);
 xtmp=NaN(NrNo,1);
 xtmp(1)=norm(model.Nodes(2,2:3)-model.Nodes(1,2:3));
 xtmp(NrNo)=norm(model.Nodes(end,2:3)-model.Nodes(end-1,2:3));
 xtmpv=model.Nodes(3:end,2:3)-model.Nodes(1:end-2,2:3);
 xtmp(2:end-1)=sqrt(sum(xtmpv.*xtmpv, 2));
 wNodes=xtmp/sum(xtmp);
 
 %% solve EigvalueProblem

 
 
 if usejava('jvm'); waitbar(0,wbrEP,'runEigenProblem EigvalueProblem');end
 lasti=length(matches);
 DisplAbs=NaN(lasti,1);
 d2ksi=NaN(1,lasti);
 if ~displacementsenable
  lasti=0;
 end
 for i = 1:lasti %i=0 keine last; i=1 erste Last
  if usejava('jvm'); waitbar(i/length(matches),wbrEP,'runEigenProblem EigvalueProblem');end
  %disp('Lambda:');
  %disp(lambda(matches(i)));
  %Lambda=fulllambda(matches(i)) %#ok<NASGU,NOPRT>
  %disp('-------------');
  %     Energy(i,2) = membrane(matches(i));
  %     Energy(i,3) = nonmembrane(matches(i));
  %for i
  if i==0
   
   
  else  % if i~= 1
   
   
   % -1.
   
   if matches(i)>=1 && displacementsenable
    if matches(i)==1
     DisplAbs(i)=mean(sqrt(sum(Displ{matches(i)}.*Displ{matches(i)},1))); % Mittelwert von [sqrt(x²+y²+z²) für jeden Knoten]
     DisplMatip1=Displ{matches(i+1)};
     DisplAbs(i+1)=mean(sqrt(sum(DisplMatip1.*DisplMatip1,1))); % Mittelwert von [sqrt(x²+y²+z²) für jeden Knoten]
     dksiMat=Displ{matches(i)+1}-Displ{matches(i)};
     dl=(model.lambda(matches(i)+1)-model.lambda(matches(i)));
     d2ksi(i)=NaN;
    elseif matches(i)==matches(end)
     dksiMat=(Displ{matches(i)}-Displ{matches(i)-1});
     dl=(model.lambda(matches(i))-model.lambda(matches(i)-1));
     d2ksi(i)=NaN;
    elseif displacementsenable
     DisplMatip1=Displ{matches(i+1)};
     DisplAbs(i+1)=mean(sqrt(sum(DisplMatip1.*DisplMatip1,1))); % Mittelwert von [sqrt(x²+y²+z²) für jeden Knoten]
     dksiMat=(Displ{matches(i)}-Displ{matches(i)-1})/2;
     dl=(model.lambda(matches(i)+1)-model.lambda(matches(i)-1))/2;
     d2ksi(i)=(DisplAbs(i+1)-2*DisplAbs(i)+DisplAbs(i-1))/(dl*dl);
    end
    dksi11vec = sqrt(sum(dksiMat.*dksiMat,1)); % sqrt(x²+y²+z²) für jeden Knoten
    if numel(dksi11vec)==numel(wNodes)
     %
    else
     warning('MyProgram:Disp','the Displ have wrong size')
     wNodes=ones(numel(dksi11vec),1)/numel(dksi11vec);
    end
    dksi11 = dot(dksi11vec,wNodes);
    displacements_ = Displ{matches(i)+1};
    %displacements_(abs(displacements_)<=1e-31)=0;%remove numeric issues close to zero
    %displacements_(abs(displacements_)<=1e-14)=0;%remove numeric issues close to zero
    dxidli=dksi11/dl;
    dUdli=max(abs(dksiMat(:)))/dl;
    dwdli=max(abs(dksiMat(end,:)))/dl;
    dvdli=max(abs(dksiMat(2,:)))/dl;
    dudli=max(abs(dksiMat(1,:)))/dl;
    duEdli=(max((dksiMat(1,:)))-min((dksiMat(1,:))))/dl;
   end
   
   % 2.
   if displacementsenable
    dksiMat=Displ{matches(i)+1}-Displ{matches(i)};
    dksi12 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
   else
    dksi12 = NaN;
   end
  end
  
  darclengths{i} = [dksi11,dksi12];
  dxidl{i}=dxidli;
  dUdl{i}=dUdli;
  dwdl{i}=dwdli;
  dvdl{i}=dvdli;
  dudl{i}=dudli;
  duEdl(i)=duEdli;
  
  displacements{i} = displacements_;
  lengthMaxJK{i}=max(abs(displacements_(:)));
  wMaxJK{i}=max(abs(displacements_(end,:)));
  vMaxJK{i+1}=max(abs(displacements_(2,:)));
  uMaxJK{i+1}=max(abs(displacements_(1,:)));
  [s1,idx1(i+1)]=max(displacements_(1,:));
  [s2,idx2(i+1)]=min(displacements_(1,:));
  sDiffJK(i+1)=s1-s2;
  arclengthHM{i}=sqrt(sum(sum(displacements_.*displacements_)));
   
   if strcmp(model.filename(1:4),'ecc-')
    if strcmp(model.filename(1:7),'ecc-B32')
     tmp=size(displacements_,2);
     LastNode=(tmp+1)/2;
     %del=[1,LastNode,LastNode+1,tmp];
     %keep=[2:LastNode-1,LastNode+2:tmp-1];
     EndDisp=displacements_(:,2:LastNode-1);
     MidDisp=displacements_(:,LastNode+2:tmp-1);
     BeamDisp=[EndDisp,MidDisp];
    else
     BeamDisp=displacements_(:,2:end-1);
    end
   else
    BeamDisp=displacements_;
   end
   Verschiebungen=sqrt(sum(BeamDisp.*BeamDisp,1));
   tmpNaNVer=isnan(Verschiebungen);
   if any(tmpNaNVer) && ~all(tmpNaNVer)
     Verschiebungen(tmpNaNVer)=0;
   end
   
   if strcmp(model.filename(1:3),'ecc')
    if i==1
     LagerUnten=model.BCMatlab(1,3);
     LagerOben=model.BCMatlab(end,3);
     KnotenGes=size(displacements_,2);
    end
    %fa=Verschiebungen(1:LastNode-1);
    %fb=Verschiebungen(2:LastNode);
    %fm=Verschiebungen(LastNode+1:end);
    ga=1;
    gb=1;%zusätliches Gewicht für Mittelknoten
    if strcmp(model.filename(1:5),'ecc-B')
    if i==1 %due to reduction
     LagerUnten=model.BCMatlab(1,3)-1;
     LagerOben=model.BCMatlab(end,3)-1;
     KnotenGes=size(displacements_,2)-2;
    end
     if strcmp(modelprops.elementtype(1:3),'B32') && modelprops.ecc~=0
      % ul UA (x-1)*S OA ol; zul x*ZS zol ->reduziert zu-> UA (x-1)*S OA
      tmp=size(BeamDisp,2);
      OberesAuflager=(tmp+1)/2;
      assert(LagerOben==OberesAuflager,'Auflager nicht erkannt')
     end
     %BeamDisp(:,1)
     %BeamDisp(:,OberesAuflager)
     arclengthJK{i+1}=(ga*sum(Verschiebungen)+gb*sum(Verschiebungen([2:LagerOben-1,LagerOben+1:end])))/((ga+gb)*size(Verschiebungen,2)-2*gb);
%      if i+1==221
%       disp(arclengthJK{i+1});
%      end
     %arclengthJK{i+1}=mean(Verschiebungen);
    elseif  strcmp(model.filename(1:7),'ecc64-B')
     if i==1
      KnotenMitte=LagerOben-64;
      if strcmp(modelprops.elementtype(1:3),'B32')
       % ul UA (x-1)*US M 63*OS OA ol; zul x*ZUS 64*ZOS zol
       StabUnten=[LagerUnten:KnotenMitte,LagerOben+3:KnotenGes-65];
       StabOben=[KnotenMitte:LagerOben,KnotenGes-64:KnotenGes-1];
      else
       StabUnten=LagerUnten:KnotenMitte;
       StabOben=KnotenMitte:LagerOben;
      end
     end
     arcUnten=(sum(Verschiebungen(StabUnten))-sum(Verschiebungen([KnotenMitte,LagerUnten]))/2)/(numel(StabUnten)-1);
     arcOben=(sum(Verschiebungen(StabOben))-sum(Verschiebungen([KnotenMitte,LagerOben]))/2)/(numel(StabOben)-1);
     arclengthJK{i+1}=(arcOben+arcUnten)/2;
%      if i+1==221
%       disp(arclengthJK{i+1});
%      end
     %arclengthJK{i+1}=mean(Verschiebungen);
    end
   else
    arclengthJK{i+1}=(sum(Verschiebungen)+2*sum(Verschiebungen(2:end-1)))/(3*size(Verschiebungen,2)-4);
    %arclengthJK{i+1}=mean(Verschiebungen);
   end
   
   



 end%for i = 1:length(matches)
 arclengthJK{1}=0*arclengthJK{2};
 vMaxJK{1}=0*vMaxJK{2};
 uMaxJK{1}=0*uMaxJK{2};
 dxidl{1}=NaN*dxidl{2};
 
 model.darclengthsJK = darclengths;
 model.arclengthuJK = arclengthJK;
 model.lengthMaxJK = lengthMaxJK;
 model.wMaxJK = wMaxJK;
 model.vMaxJK = vMaxJK;
 model.uMaxJK = uMaxJK;
 model.dUdl=dUdl;
 model.dwdl=dwdl;
 model.dvdl=dvdl;
 model.dudl=dudl;
 model.duEdl=duEdl; 
 model.arclengthuHM = arclengthHM;
 model.displacementsJK = displacements;
 model.dxidl=dxidl;
 model.d2ksi=d2ksi;
 model.uDiffJK = sDiffJK;%veraltet use model.uDiffJK 
 model.sDiffJK = sDiffJK;
 model.idx1 = idx1;
 model.idx2 = idx2;
 

 
end %fucntion

