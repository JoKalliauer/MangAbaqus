function model = runEigenProblemDispJK(~,model,Displ,~,~,matches,wbrEP)
 
 % [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);
 
 %Displ = NodalResults2Displ(Nres);
 if numel(Displ)<1
  displacementsenable=false;
 else
  displacementsenable=true;
 end
 %Kg = EigRes;
 lambda0=model.lambda0';

 displacements = cell(length(lambda0),1);
 darclengths = cell(length(lambda0),1);
 arclengthJK = cell(length(lambda0),1);
 lengthMaxJK = cell(length(lambda0),1);
 wMaxJK = cell(length(lambda0),1);
 vMaxJK = cell(length(lambda0),1);
 arclengthHM = cell(length(lambda0),1);
 dxidl = cell(length(lambda0),1);
 dUdl = cell(length(lambda0),1);
 dwdl = cell(length(lambda0),1);
 dvdl = cell(length(lambda0),1);
 
 NrNo=size(model.Nodes,1);
 xtmp=NaN(NrNo,1);
 xtmp(1)=norm(model.Nodes(2,2:3)-model.Nodes(1,2:3));
 xtmp(NrNo)=norm(model.Nodes(end,2:3)-model.Nodes(end-1,2:3));
 xtmpv=model.Nodes(3:end,2:3)-model.Nodes(1:end-2,2:3);
 xtmp(2:end-1)=sqrt(sum(xtmpv.*xtmpv, 2));
 wNodes=xtmp/sum(xtmp);
 
 %% solve EigvalueProblem

 
 
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


   % -1.
   if matches(i)>1
    dksiMat=Displ{matches(i)}-Displ{matches(i)-1};
    dksi11vec = sqrt(sum(dksiMat.*dksiMat,1));
    if numel(dksi11vec)==numel(wNodes)
     %
    else
     warning('MyProgram:Disp','the Displ have wrong size')
     wNodes=ones(numel(dksi11vec),1)/numel(dksi11vec);
    end
    dksi11 = dot(dksi11vec,wNodes);
    displacements_ = Displ{matches(i)-1};
    %displacements_(abs(displacements_)<=1e-31)=0;%remove numeric issues close to zero
    %displacements_(abs(displacements_)<=1e-14)=0;%remove numeric issues close to zero
    dl=(model.lambda(matches(i))-model.lambda(matches(i)-1));
    dxidli=dksi11/dl;
    dUdli=max(abs(dksiMat(:)))/dl;
    dwdli=max(abs(dksiMat(end,:)))/dl;
    dvdli=max(abs(dksiMat(2,:)))/dl;
   else
    dksi11 = NaN;
    displacements_(:,:) = NaN(3,1);
    dxidli=NaN;
    dUdli=NaN;
    dwdli=NaN;
    dvdli=NaN;
   end

   % 2.
   if displacementsenable
    dksiMat=Displ{matches(i)+1}-Displ{matches(i)};
    dksi12 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
   else
    dksi12 = NaN;
   end

   % 4.
%    if displacementsenable
% %    displacements_(:,:,6) = Displ{matches(i)-0};
% %    displacements_(:,:,7) = Displ{matches(i)+1};
% %     displacements_(:,:,8) = Displ{matches(i)+2};
%    end
   



   end
   
%    if strcmp(modelprops.typeofanalysis,'KNL2') || mintest<matches(i)+4
% %     if displacementsenable; displacements_(:,:,9)=NaN*displacements_(:,:,8);end
%    else
%     %displacements_(:,:,9) = Displ{matches(i)+3};
%    end
   
   darclengths{i} = [dksi11,dksi12];
   dxidl{i}=dxidli;
   dUdl{i}=dUdli;
   dwdl{i}=dwdli;
   dvdl{i}=dvdli;
   
   displacements{i} = displacements_;
   arclengthJK{i}=mean(sqrt(sum(displacements_.*displacements_,1)));
   lengthMaxJK{i}=max(abs(displacements_(:)));
   wMaxJK{i}=max(abs(displacements_(end,:)));
   vMaxJK{i}=max(abs(displacements_(2,:)));
   arclengthHM{i}=sqrt(sum(sum(displacements_.*displacements_)));


 end%for i = 1:length(matches)
 
 model.darclengthsJK = darclengths;
 model.arclengthuJK = arclengthJK;
 model.lengthMaxJK = lengthMaxJK;
 model.wMaxJK = wMaxJK;
 model.vMaxJK = vMaxJK;
 model.dUdl=dUdl;
 model.dwdl=dwdl;
 model.dvdl=dvdl;
 model.arclengthuHM = arclengthHM;
 model.displacementsJK = displacements;
 model.dxidl=dxidl;

 

 
end %fucntion

