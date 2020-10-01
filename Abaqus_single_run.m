%#!/bin/rm
%university:TU Wien

function [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main)
 %% Inputkontrolle
 assert(modelprops.numofelm>0,'numberOfElements must be greater zero')
 mustBeInteger(modelprops.numofelm)
 modelprops.numofelm=uint16(modelprops.numofelm);
 assert(numel(modelprops.lambda)<41,'using %f>27 lambda-values takes much resources',numel(modelprops.lambda))
 modelprops.numofeigs=min(modelprops.numofeigs,uint8(14));
 %%% Setting default values

 
 %% Programm
 model = runEigenProblem(modelprops);
 model.check=main.check;
 model.filename=[model.filename,'-',modelprops.typeofanalysis];
 for i=1:numel(model.eigenvalues)
  if numel(model.eigenvalues{i})==0
   warning('MyProgram:Empty','some eigenvalues are empty, please rerun using modelprops.forceAbaqus=true')
   model.eigenvalues=model.eigenvalues(1:i-1);
   model.lambda0=model.lambda0(1:i-1);
   break
  end
  if ~isreal(model.eigenvalues{i}(5:6,:))
   warning('MyProgram:Komplex','Some values are komplex %f',i)
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
 end
 %JKExtend(model,modelprops);
 if sum(strcmp(fieldnames(modelprops), 'savefigures')) == 1
  model.savefigures=modelprops.savefigures;
 else
  model.savefigures=false;
 end
 if modelprops.closall==true
  close all
 end
 if (strcmpi(sortType,'none'))&&isempty(forcedeig)
  NrEWs=min([modelprops.numofeigs,size(model.eigenvalues{1},2)]);
  if NrEWs>7
   warning('Myprogram:color','Be aware plotting more than 7 graphs might lead to the same/repeating colour')
  end
  for k3 = 1:1:NrEWs
   res = sortEigenValuesAndGetQuantities(model,sortType,plotfig,k3);
  end
 else
  res = sortEigenValuesAndGetQuantities(model,sortType,plotfig,forcedeig);
 end
end
