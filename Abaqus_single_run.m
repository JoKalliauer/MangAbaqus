function [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,numofelm,ecc,elementtype)
%% Thats the function to call, which calls all the required subfunctions
%Outermost function that calls all programs, but Input must be provieded by a skript
% modelprops.* enter runEigenProblem, hoever main.* does not enter runEigenProblem and is for postprocess, model.* is created by Matlabfunctions, not by the user.
% e.g. [model] = runEigenProblem(modelprops)
%university:TU Wien
%author: Johannes Kalliauer, Michal Malendowski

%% Input
% modelpros.* Input for runEigenProblem
%  modelpros.numofelm numer of elements, might get overwritten by "numofelm"
% sortType ... How the eigenvalues/-vectors are sorted
% plotfig  ... which figues to plot
% forcedeig... Evaluate thos eigenvalues
% main.* Inputvariabels for postprocess
% numofelem ... Number of Elements
% ecc       ... Eccentricity
% elementtype...Type of element

%% Output
% res  ... resulst of postprocess
% model... results of runEigenProblem

%% Structure
% * runEigenProblem ... run the Eingenvalue-Problem
%   * selectModel .. calls a function to create the input-file 
%   * AbaqusModelsGeneration.runAbaqus ... run the input-file in Abaqus
%   * AbaqusModelsGeneration.getStiffnessMatrices ... get the stiffness-matrix from Abaqus-results
%   * AbaqusModelsGeneration.getHistoryOutputFromDatFile ... get the nodal-results from Abaus
%   * runEigenProblemSub ... Run the core of the eigenvalue-Problem and saving it into model
%     * solveCLEforMinEigNew ... get one specific eigenvalue and the eigenvector
%   * runEigenProblemDispJK ... Posprocessing the displacements
% * sortEigenValuesAndGetQuantities ... does the calculation of \rho
% * plotresMulti ... Plots the requested graphs

%% recent-change
%2023-02-16 JK: enabling sound

%% Code

warning('on','MyProgramm:lowPrecission:RhoOne')

%% Inputkontrolle
if exist('ecc','var')
 modelprops.ecc=ecc;
end
if sum(strcmp(fieldnames(modelprops), 'numelFac')) == 0
 modelprops.numelFac=1;
 numelFac=1;
else
 numelFac=modelprops.numelFac;
end
if exist('numofelm','var')
 if numel(numofelm)>0
  modelprops.numofelm=numofelm*numelFac;
 end
end
if exist('elementtype','var')
 modelprops.elementtype=elementtype;
end
assert(min(modelprops.numofelm)>0,'numberOfElements must be greater zero')
mustBeInteger(modelprops.numofelm)
modelprops.numofelm=uint16(modelprops.numofelm);
if sum(strcmp(fieldnames(modelprops), 'numofeigs')) == 0
 modelprops.numofeigs=1; %setting default value to 1
end
if modelprops.numofeigs>56%14
 warning('MyProgramm:Input','For higher precission reduce the number of requested eigenvalues below 14')
 if modelprops.numofeigs>2500
  modelprops.numofeigs=min(2500,20*(modelprops.numofelm+2));
  warning('MyProgramm:Input','The eigenvalues are reduced to %d, otherwise it will need more than 30GB RAM',modelprops.numofeigs)
 end
end
%  if any(modelprops.lambda<4*modelprops.epsilon)
%   warning('MyProgram:Input','first lamdavalue must be four times epsilon or larger')
%   modelprops.lambda=sort(modelprops.lambda(modelprops.lambda>=3*modelprops.epsilon));
%  end
%assert(all(modelprops.lambda>=4*modelprops.epsilon),'first lamdavalue must be four times epsilon or larger')
if numel(modelprops.lambda)>301
 warning('MyProgram:Input','using %f>201 lambda-values takes much resources',numel(modelprops.lambda))
 if numel(modelprops.lambda)>=331
  warning('MyProgram:Input','using %f>=331 Abaqus CAE does not load the stucture',numel(modelprops.lambda)) %ecc-B32H-20-l5-e0.040447-f1-eps0.02-u1
 end
 %assert(numel(modelprops.lambda)<=1001,'using %f>401 lambda-values takes much resources',numel(modelprops.lambda))
 %  assert(numel(modelprops.lambda)<=2195,'using %f>401 lambda-values takes much resources',numel(modelprops.lambda))
end
assert(max([forcedeig 0])<=min(modelprops.numofeigs),'forceeig must be smaler or equal that number of calculated ones')
modelprops.forcedeig=forcedeig;
[sl1, sl2]=size(modelprops.lambda);
if sl1>1 && sl2==1; modelprops.lambda=transpose(modelprops.lambda); end
modelprops.lambda=unique(sort([modelprops.lambda]));
diffs=modelprops.lambda(2:end) - modelprops.lambda(1:end-1);
mindiff=min(diffs);
if numel(mindiff)>0
 assert(modelprops.epsilon-mindiff<=+eps(2*modelprops.lambda(end)),'epsilon larger than lamdba0-steps')
end
if mindiff<7*modelprops.epsilon
 maxdiff=max(diffs(2:end));
 if mindiff<maxdiff+eps(max(1,maxdiff))
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
if sum(strcmp(fieldnames(main), 'rstabil')) == 0
 main.rstabil=0.9999;
end
if modelprops.numofelm>1024
%assert(modelprops.numofelm<=1024,'more than 1000 Elements need more than 24GB RAM')
%assert(sum(modelprops.numofelm)<=2000,'more than 1000 Elements need more than 24GB RAM')
%assert(sum(modelprops.numofelm)<=2048,'more than 1000 Elements need more than 24GB RAM')
assert(sum(modelprops.numofelm)<=2*2048,'more than 1000 Elements need more than 24GB RAM')
 if strcmp(modelprops.elementtype(end),'H')
  assert(sum(modelprops.numofelm)<2000,'MATLAB-Problem with eigs: Maximum number of attempts to perform numeric factorization of symmetric indefinite matrix exceeded.')%BB5-B31OSH2000-l5-f1-eps0.001-u1
 end
end
if modelprops.numofelm*numel(modelprops.lambda)>=131070
 assert(modelprops.numofelm<200,'Matlab will need more than 32GB')
 assert(numel(modelprops.lambda)<671,'Matlab will need more than 32GB')
end
if sum(strcmp(fieldnames(main), 'whichEV')) == 0
 main.whichEV='split';
end
modelprops.whichEV=main.whichEV;
if sum(strcmp(fieldnames(main), 'Normierung')) == 0
 main.Normierung='rCT_K0_r';
end
modelprops.Normierung=main.Normierung;
if strcmp(main.whichEV,'split')
 main.whichEV={'Disp'};
elseif modelprops.numofeigs==0
 main.whichEV='split';
end
if strcmp(main.whichEV,'Hyb') && modelprops.numofeigs>0
 assert(strcmp(modelprops.elementtype(end),'H'),'modelprops.elementtype does not have hybrid dofs')
end
modelprops.whichEV=main.whichEV;
if ismember(943,plotfig)
 if strcmp(modelprops.whichEV,'bungle_rKr') || strcmp(modelprops.whichEV,'bungle_rK0r') || strcmp(modelprops.whichEV,'bungle')
  %
 else
  assert(0,'plot943 is proposed to work with bungle_rKr')
 end
end
if sum(strcmp(fieldnames(main), 'savefigures')) == 0
 main.savefigures=false;
end
if sum(strcmp(fieldnames(main), 'xBezug')) == 0
 main.xBezug='n';
end
if sum(strcmp(fieldnames(main), 'flipAxis')) == 0
 main.flipAxis=false;
end



%%% Setting default values


%% Programm

if modelprops.loadfactor==0 %&& numel(modelprops.numelFac)>1
 modelpropsNEW=modelprops;
 modelpropsNEW.loadfactor=1;
 [modelP1] = runEigenProblem(modelpropsNEW);
 if numel(modelprops.numelFac)>1
  %modelpropsNEW.numelFac=fliplr(modelpropsNEW.numelFac);
  modelpropsNEW.numofelm=fliplr(modelpropsNEW.numofelm);
 else
  modelpropsNEW.loadfactor=-1;
 end
 [modelM1] = runEigenProblem(modelpropsNEW);
 modelpropsNEW.loadfactor=NaN;
 modelpropsNEW.numofelm=NaN;
 model=mergeModel(modelP1,modelM1);
else
 %%%%%%%%%%%%%%  [model] = runEigenProblem(modelprops)  %%%%%%%%%%%%%%
 [model] = runEigenProblem(modelprops);
end
%model.fullEV

%% post process
%  main.numofeigs=modelprops.numofeigs;
model.check=main.check;
model.filename=strcat(model.filename,'-',modelprops.typeofanalysis,'-',char(main.whichEV));
model.lambdainput=model.lambda0;
main.allowComplex=modelprops.allowComplex;
for i=1:numel(model.eigenvalues)
 if numel(model.eigenvalues{i})==0
  if modelprops.numofeigs>0 && i<numel(model.eigenvalues)
   warning('MyProgram:Empty','some eigenvalues are empty, please rerun using modelprops.forceAbaqus=true')
  end
  model.eigenvalues=model.eigenvalues(1:i-1);%delete last entry
  %model.lambda0=model.lambda0(1:i-1);
  break
 end
 if ~isreal(model.eigenvalues{i}(5,:)) && ~strcmp(modelprops.whichEV,'skip')
  %warning('MyProgram:Komplex','Some values are komplex %f',i)
  %model.eigenvalues{i}=real(model.eigenvalues{i});
  EW=model.eigenvalues{i}(5,:);
  wo=(abs(EW)-abs(real(EW)))>4;
  if any(wo)
   [LaststufenJK,DOFsJK,~]=size(model.eigenvectors{i});
   model.eigenvectors{i}(:,:,wo)=NaN(LaststufenJK,DOFsJK,sum(wo)); %(Last,DOF,EW)
  else
   if sum(strcmp(fieldnames(model), 'eigenvectors')) == 0
    error('MyPrgm:Rerun','Try reruning with forcerun=1 and maybe with different whichEV')
   end
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
if numel(model.eigenvalues)>=1
 NrEWs=min([modelprops.numofeigs,size(model.eigenvalues{1},2)]);
else
 NrEWs=modelprops.numofeigs;
end
if isempty(forcedeig)
 resEWs=1:NrEWs;
 if NrEWs>7
  if NrEWs>19
   warning('Myprogram:color','Be aware plotting more than 19 graphs will definitly lead to the repeating colours, you are using %d',NrEWs)
  else
   %warning('Myprogram:color','Be aware plotting more than 7 graphs might lead to the same/repeating colour')
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



lastLambda=max(modelprops.lambda);
minLambda=min([modelprops.lambda,-lastLambda]);
main.typeofanalysis=modelprops.typeofanalysis;
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
 if strcmp(modelprops.whichEV,'skip')
  res=NaN;
 else
  %%%%%%%%%%%%%%  sortEigenValuesAndGetQuantities  %%%%%%%%%%%%%%
  res(k3) = sortEigenValuesAndGetQuantities(model,sortType,[],k3,limit,lastLambda,main);  % forcedeig=k3;
  toolarge=res(k3).lambda>max(modelprops.lambda);
  if any(toolarge)
   warning('MyProgram:Input:LamdbaAvail:SingleRun','more lamdas available than requested, extend the lambda-range or try using modelprops.forcerun=true')
   warning('off','MyProgram:Input:LamdbaAvail:SingleRun')
   res(k3).lambda(toolarge)=NaN;
  end
  toostart = abs(res(k3).lambda(2:end)) < (min(modelprops.lambda(1:end))-eps(single(10)));
  if any(toostart)
   warning('MyProgram:Input','first lambdas ignored try using modelprops.forcerun=true')
   res(k3).lambda(toostart)=NaN;
  end
 end
end % for k3 = resEWs
if numel(resEWs)==0
 res=NaN;
end
toolarge = abs(model.lambda0)>max(abs(modelprops.lambda))+eps(single(1));
if any(toolarge)
 warning('MyProgram:Input','more lamdas available than requested try using modelprops.forcerun=true')
 model.lambda0(toolarge)=NaN;
 model.DetKtx(toolarge)=NaN;
 model.load0(toolarge)=NaN;
end
toolarge=model.lambda>max(modelprops.lambda);
if any(toolarge)
 warning('MyProgram:Input','more lamdas available than requested try using modelprops.forcerun=true')
 model.lambda(toolarge)=NaN;
 %   model.DetKtx(toolarge)=NaN;
 model.load(toolarge)=NaN;
 model.load0([false;toolarge])=NaN;
end
toolarge=model.lambda0>max(modelprops.lambda)+eps(single(1));
if any(toolarge)
 warning('MyProgram:Input','more lamdas available than requested try using modelprops.forcerun=true')
 model.lambda0(toolarge)=NaN;
 %   model.DetKtx(toolarge)=NaN;
 model.load0(toolarge)=NaN;
end
toolarge=model.fulllambda>max(modelprops.lambda);
if any(toolarge)
 model.fulllambda(toolarge)=NaN;
 model.fullload0(toolarge)=NaN;
end
toonegativ=model.fulllambda<-max(modelprops.lambda);
if any(toonegativ)
 model.fulllambda(toonegativ)=NaN;
 model.fullload0(toonegativ)=NaN;
end
toostart=model.fulllambda<min(modelprops.lambda);
if any(toostart)
 model.fulllambda(toostart)=NaN;
end
toostart=model.lambda<min(modelprops.lambda);
if any(toostart)
 model.lambda(toostart)=NaN;
end


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
 if strcmp(main.whichEV,'skip') && size(model.eigenvalues,1)~=0
  lengthAbaqus=size(model.eigenvalues,1)+4;
 else
  lengthAbaqus=size(model.fullEV,2);
 end
 lengthInput=size(model.fulllambda,1);
 if lengthInput>lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forcerun=true')
  model.fulllambda=model.fulllambda(1:lengthAbaqus);
 elseif lengthInput<lengthAbaqus
  warning('MyProgram:Input','lamdalengh differs, try modelprops.forceAbaqus=true')
  model.fullEV=model.fullEV(:,1:lengthInput);
 end
%  if any(all(isnan(real(model.fullEV)))) %wenn es Spalten gibt die in jeder Zeile immer NaN sind
%   model.fulllambda(all(isnan(real(model.fullEV)), 1)) = [];
%   %model.fullEV(:,all(isnan(real(model.fullEV)), 1)) = [];
%  end
 tooLarge= (model.fulllambda>lastLambda+eps(lastLambda));
 if any(tooLarge)
  model.fulllambda(tooLarge)=NaN;
  tooNegativ=(model.fulllambda<minLambda-eps(minLambda));
  model.fulllambda(tooNegativ)=NaN;
 end
 %%%%%%%%%%%%%%  plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main)  %%%%%%%%%%%%%%
 plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main);
end%if usejava('desktop')
if main.savefigures==true %&& ~strcmp(main.whichEV,'skip')
 printresMulti(res,model,plotfig,[],[],resEWs,main.whichEV)
end
disp(['finish: ','AnalysisResults/',model.filename,'-',num2str(model.numofeigs),'.mat']);

  %sound when finished (easter egg)
  beep
  load handel y
  sound(y,8192) %remove this line if Device Error: Illegal combination of I/O devices

end
