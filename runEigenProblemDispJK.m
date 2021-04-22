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
 arclengthHM = cell(length(lambda0),1);
 dxidl = cell(length(lambda0),1);
 
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
   % -4.
   Zeile=matches(i)-4;
   if Zeile>0
    if Zeile>1
     dksiMat=Displ{Zeile}-Displ{matches(i)-5};
     dksi04 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
%      displacements_(:,:,1) = Displ{matches(i)-5};
    else
     dksi04 = NaN;
%      displacements_(:,:,1) = NaN*Displ{1};
    end
    dksiMat=Displ{matches(i)-3}-Displ{matches(i)-4};
    dksi03 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
%     displacements_(:,:,2) = Displ{matches(i)-4};
   else
    %Kt04 = 0*Kts{1,2};
    dksi04 = NaN;
    dksi03 = NaN;
%     if displacementsenable
%      displacements_(:,:,1) = NaN*Displ{1};
%      displacements_(:,:,2) = NaN*Displ{1};
%     end
   end
   %Kt04(ru,:) = []; Kt04(:,ru) = [];
   % -3.
   if Zeile>-1
    dksiMat=Displ{matches(i)-2}-Displ{matches(i)-3};
    dksi02 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
%     displacements_(:,:,3) = Displ{matches(i)-3};
   else
    dksi02 = NaN;
%     displacements_(:,:,3) = NaN;
   end
   % -2.
   if Zeile>-2
    dksiMat=Displ{matches(i)-1}-Displ{matches(i)-2};
%     dksi01vec = sqrt(sum(dksiMat.*dksiMat,1));
    dksi01 = mean(sqrt(sum(dksiMat.*dksiMat,1)));
%     displacements_(:,:,4) = Displ{matches(i)-2};
   else
    dksi01 = NaN;
%     displacements_(:,:,4) = NaN;
   end

   % -1.
   if Zeile>-3
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
    dxidli=dksi11/(model.lambda(matches(i))-model.lambda(matches(i)-1));
   else
    dksi11 = NaN;
    displacements_(:,:) = NaN;
    dxidli=NaN;
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
   
   dksi = [dksi04, dksi03, dksi02, dksi01, dksi11,dksi12];
   darclengths{i} = dksi;
   dxidl{i}=dxidli;
   
   displacements{i} = displacements_;
   arclengthJK{i}=mean(sqrt(sum(displacements_.*displacements_,1)));
   arclengthHM{i}=sqrt(sum(sum(displacements_.*displacements_)));


 end%for i = 1:length(matches)
 
 model.darclengthsJK = darclengths;
 model.arclengthuJK = arclengthJK;
 model.arclengthuHM = arclengthHM;
 model.displacementsJK = displacements;
 model.dxidl=dxidl;

 

 
end %fucntion

