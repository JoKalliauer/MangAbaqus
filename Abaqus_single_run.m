%#!/bin/rm
%university:TU Wien

<<<<<<< HEAD
function [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,ecc)

%% Inputkontrolle
 if exist('ecc','var')
  modelprops.ecc=ecc;
 end
 assert(min(modelprops.numofelm)>0,'numberOfElements must be greater zero')
=======
function [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main)

%% Inputkontrolle
 assert(modelprops.numofelm>0,'numberOfElements must be greater zero')
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 mustBeInteger(modelprops.numofelm)
 modelprops.numofelm=uint16(modelprops.numofelm);
 %modelprops.numofeigs=min(modelprops.numofeigs,uint8(14));
 if modelprops.numofeigs>14
  warning('MyProgramm:Input','For higher precission reduce the number of requested eigenvalues')
 end
%  if any(modelprops.lambda<4*modelprops.epsilon)
%   warning('MyProgram:Input','first lamdavalue must be four times epsilon or larger')
%   modelprops.lambda=sort(modelprops.lambda(modelprops.lambda>=3*modelprops.epsilon));
%  end
 %assert(all(modelprops.lambda>=4*modelprops.epsilon),'first lamdavalue must be four times epsilon or larger')
 if numel(modelprops.lambda)>201
  warning('MyProgram:Input','using %f>201 lambda-values takes much resources',numel(modelprops.lambda))
  %assert(numel(modelprops.lambda)<=1001,'using %f>401 lambda-values takes much resources',numel(modelprops.lambda))
%  assert(numel(modelprops.lambda)<=2195,'using %f>401 lambda-values takes much resources',numel(modelprops.lambda))
 end
<<<<<<< HEAD
 assert(max([forcedeig 0])<=modelprops.numofeigs,'forceeig must be smaler or equal that number of calculated ones')
 modelprops.lambda=unique(sort([0 modelprops.lambda]));
 diffs=modelprops.lambda(2:end) - modelprops.lambda(1:end-1);
 mindiff=min(diffs);
 if numel(mindiff)>0
  assert(modelprops.epsilon-mindiff<=+eps(2*modelprops.lambda(end)),'epsilon larger than lamdba0-steps')
 end
 if mindiff<7*modelprops.epsilon
  maxdiff=max(diffs);
  if mindiff<maxdiff+eps(max(1,maxdiff))
=======
 assert(max([forcedeig 1])<=modelprops.numofeigs,'forceeig must be smaler or equal that number of calculated ones')
 modelprops.lambda=unique(sort([0 modelprops.lambda]));
 diffs=modelprops.lambda(2:end) - modelprops.lambda(1:end-1);
 mindiff=min(diffs);
 assert(modelprops.epsilon-mindiff<=+eps(2*modelprops.lambda(end)),'epsilon larger than lamdba0-steps')
 if mindiff<7*modelprops.epsilon
  maxdiff=max(diffs);
  if mindiff<maxdiff+eps(1)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
   ratio=maxdiff/modelprops.epsilon;
   assert(abs(ratio-uint8(ratio))<eps(35),'epsilon larger than lamdba0-steps')
  else
   warning('MProgram:Input','Probably overlapping full-lamda-values')
  end
 end
 plotfig=unique(sort(plotfig));% just to be easier to read for humans
 if sum(strcmp(fieldnames(modelprops), 'forceAbaqus')) == 0
  modelprops.forceAbaqus=false;
 elseif modelprops.forceAbaqus==true
  %modelprops.forcerun=true;
 end
<<<<<<< HEAD
 if sum(strcmp(fieldnames(main), 'rstabil')) == 0
   main.rstabil=0.9999;
 end
 %assert(modelprops.numofelm<=1024,'more than 1000 Elements need more than 24GB RAM')
 assert(sum(modelprops.numofelm)<=2000,'more than 1000 Elements need more than 24GB RAM')
=======
 assert(modelprops.numofelm<=1024,'more than 1000 Elements need more than 24GB RAM')
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
 %%% Setting default values

 
 %% Programm
 
 model = runEigenProblem(modelprops);
 %model.fullEV
 
 %% post process

 model.check=main.check;
 model.filename=[model.filename,'-',modelprops.typeofanalysis];
 model.lambdainput=model.lambda0;
 for i=1:numel(model.eigenvalues)
  if numel(model.eigenvalues{i})==0
   warning('MyProgram:Empty','some eigenvalues are empty, please rerun using modelprops.forceAbaqus=true')
   model.eigenvalues=model.eigenvalues(1:i-1);
   model.lambda0=model.lambda0(1:i-1);
   break
  end
<<<<<<< HEAD
  if ~isreal(model.eigenvalues{i}(5,:))
=======
  if ~isreal(model.eigenvalues{i}(5:6,:))
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
   %warning('MyProgram:Komplex','Some values are komplex %f',i)
   %model.eigenvalues{i}=real(model.eigenvalues{i});
   EW=model.eigenvalues{i}(5,:);
   wo=(abs(EW)-abs(real(EW)))>4;
   if any(wo)
    [LaststufenJK,DOFsJK,~]=size(model.eigenvectors{i});
    model.eigenvectors{i}(:,:,wo)=NaN(LaststufenJK,DOFsJK,sum(wo)); %(Last,DOF,EW)
   else
    model.eigenvectors{i}=real(model.eigenvectors{i});
   end
  end
 end %for i=1:numel(model.eigenvalues)
 %JKExtend(model,modelprops);
 if sum(strcmp(fieldnames(main), 'savefigures')) == 1
  model.savefigures=main.savefigures;
 else
  model.savefigures=false;
 end
%  model.savefigures=false;
<<<<<<< HEAD
 if numel(model.eigenvalues)>=1
  NrEWs=min([modelprops.numofeigs,size(model.eigenvalues{1},2)]);
 else
  NrEWs=modelprops.numofeigs;
 end
 if isempty(forcedeig)
=======
 NrEWs=min([modelprops.numofeigs,size(model.eigenvalues{1},2)]);
 if (strcmpi(sortType,'none'))&&isempty(forcedeig)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  resEWs=1:NrEWs;
  if NrEWs>7
   if NrEWs>19
    warning('Myprogram:color','Be aware plotting more than 19 graphs will definitly lead to the repeating colours, you are using %f',NrEWs)
   else
    warning('Myprogram:color','Be aware plotting more than 7 graphs might lead to the same/repeating colour')
   end
  end
 else
  if any(forcedeig>NrEWs)
   warning('MyProgram:Input','not enough NrEWs for some forcedeig')
   resEWs=forcedeig(forcedeig<NrEWs);
  else
   resEWs=forcedeig;
  end
 end %if (strcmpi(sortType,'none'))&&isempty(forcedeig)

<<<<<<< HEAD
 if numel(resEWs)==0
  res=NaN;
 end
 
 lastLambda=max(modelprops.lambda);
 for k3 = resEWs
  if strcmp(modelprops.testcase,'eccenCompressionBeam') && strcmp(modelprops.elementtype,'B32OS') && strcmp(modelprops.typeofanalysis,'KNL2') && modelprops.epsilon == 0.01 && max(modelprops.lambda)>1.8
   %limit.OC5a=0.49e-4;
   %limit.C5diff=-0.011;
   limit.empty=true;
  else
   limit.new=false;
  end
  if strcmp(modelprops.testcase,'pureBendingBeam') && strcmp(modelprops.elementtype,'B32OSH') && strcmp(modelprops.typeofanalysis,'KNL2')
   limit.C8minrho=.991;
  end
  res(k3) = sortEigenValuesAndGetQuantities(model,sortType,plotfig,k3,limit,lastLambda,main); 
  toolarge=res(k3).lambda>max(modelprops.lambda);
  if any(toolarge)
   warning('MyProgram:Input','more lamdas available than requested try using modelprops.forcerun=true')
   res(k3).lambda(toolarge)=NaN; 
  end
  toostart=res(k3).lambda(2:end)<min(modelprops.lambda(2:end));
  if any(toostart)
   warning('MyProgram:Input','first lambdas ignored try using modelprops.forcerun=true')
   res(k3).lambda(toostart)=NaN; 
  end
 end % for k3 = resEWs
 
 
 
if usejava('desktop')
 if main.closall==true
  if main.colorshift~=0
   warning('MyProgram:Input','Setting main.colorshift to zero')
   main.colorshift=0;
  end
  close all
 end
 
 MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
 MyMarker=['o' 'x'];
 %colJK=MyColours{mod(k3-1,19)+1};
 %markJK=MyMarker(mod(k3-1,2)+1);
 lengthAbaqus=size(model.fullEV,2);
 lengthInput=size(model.fulllambda,1);
 if lengthInput>lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forcerun=true')
  model.fulllambda=model.fulllambda(1:lengthAbaqus);
 elseif lengthInput<lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forceAbaqus=true')
  model.fullEV=model.fullEV(:,1:lengthInput);
 end
=======
 lastLambda=max(modelprops.lambda);
 for k3 = resEWs
  if strcmp(modelprops.testcase,'eccenCompressionBeam') && strcmp(modelprops.elementtype,'B32OS') && strcmp(modelprops.typeofanalysis,'KNL2') && modelprops.epsilon == 0.01 && max(modelprops.lambda)>1.8
   %limit.OC5a=0.49e-4;
   %limit.C5diff=-0.011;
   limit.empty=true;
  else
   limit.new=false;
  end
  if strcmp(modelprops.testcase,'pureBendingBeam') && strcmp(modelprops.elementtype,'B32OSH') && strcmp(modelprops.typeofanalysis,'KNL2')
   limit.C8minrho=.991;
  end
  res(k3) = sortEigenValuesAndGetQuantities(model,sortType,plotfig,k3,limit,lastLambda); %#ok<AGROW>
  toolarge=res(k3).lambda>max(modelprops.lambda);
  if any(toolarge)
   warning('MyProgram:Input','more lamdas available than requested try using modelprops.forcerun=true')
   res(k3).lambda(toolarge)=NaN; %#ok<AGROW>
  end
  toostart=res(k3).lambda<min(modelprops.lambda(2:end));
  if any(toostart)
   warning('MyProgram:Input','first lambdas ignored try using modelprops.forcerun=true')
   res(k3).lambda(toostart)=NaN; %#ok<AGROW>
  end
 end % for k3 = resEWs
 
 
 
 if usejava('desktop')
 if main.closall==true
  if main.colorshift~=0
   warning('MyProgram:Input','Setting main.colorshift to zero')
   main.colorshift=0;
  end
  close all
 end
 
 MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
 MyMarker=['o' 'x'];
 %colJK=MyColours{mod(k3-1,19)+1};
 %markJK=MyMarker(mod(k3-1,2)+1);
 lengthAbaqus=size(model.fullEV,2);
 lengthInput=size(model.fulllambda,1);
 if lengthInput>lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forcerun=true')
  model.fulllambda=model.fulllambda(1:lengthAbaqus);
 elseif lengthInput<lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forceAbaqus=true')
  model.fullEV=model.fullEV(:,1:lengthInput);
 end
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 if any(all(isnan(real(model.fullEV)))) %wenn es Spalten gibt die in jeder Zeile immer NaN sind
  model.fulllambda(all(isnan(real(model.fullEV)), 1)) = [];
  model.fullEV(:,all(isnan(real(model.fullEV)), 1)) = [];
 end
 tooLarge= (model.fulllambda>lastLambda);
 if any(tooLarge)
  model.fulllambda(tooLarge)=NaN;
 end
<<<<<<< HEAD
 main.typeofanalysis=modelprops.typeofanalysis;
 plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main)
end%if usejava('desktop')
if main.savefigures==true
 printresMulti(res,model,plotfig,[],[],resEWs)
end
 disp(['finish: ','AnalysisResults/',model.filename,'-',num2str(model.numofeigs),'.mat']);
=======
 plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main)
end
if main.savefigures==true
 printresMulti(res,model,plotfig,[],[],resEWs)
end
 disp(['finish: ','AnalysisResults/',model.filename,'-',num2str(modelprops.numofeigs),'.mat']);
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 beep
end
