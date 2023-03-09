function printresMulti(res,model,plotfig,~,~,resEWs,mainwhichEV)


  if ~exist('Output/Figures/CSV', 'dir')
   if isunix
    mkdir('Output/Figures/CSV');
   elseif ispc
    warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
    return
   elseif ismac
    warning('MyProgram:OS','Mac not tested try unix')
   end
  end
  
 if numel(resEWs)>0 && ~strcmp(mainwhichEV,'skip')
  lambda = res(min(resEWs)).lambda;
 else
  lambda=model.lambda;
 end
 lengthlam=numel(lambda);
 if lengthlam==0
  lengthlam=numel(res(resEWs(1)).lambda);
 end
 m14=NaN(lengthlam,1+numel(resEWs));
 m14(:,1)=lambda;
 m07=m14;
 %m15=m14;%actually it is fulllenght and not lengthlam
 m20=m14;
 m34=m14;
 
 loadlen=min(lengthlam-1,numel(model.load));
 m30=NaN(loadlen+1,1+numel(resEWs));
 
 
 [tmp1,tmp2]=size(model.load);
 if tmp1==1 && tmp2>1
  model.load=transpose(model.load);
 end
 m30(:,1)=[0;model.load(1:loadlen)];
 m32=m30;
 
   xPlot15=model.fulllambda(1:end);
   m15=NaN(numel(xPlot15),1+numel(resEWs));
   m15(:,1)=xPlot15;
   m19=m15;
 for k3=resEWs
  i=find(resEWs==k3);
  if ismember(7,plotfig)
   m07(:,1+i)=res(k3).TAU;
   filename=strcat('Output/Figures/CSV/',model.filename,'_tau07.csv');
   writematrix(m07,filename,'Delimiter',';')
   %system(['exec sed -i "s/\./,/g" ',filename]);
  end
  if ismember(14,plotfig) && isstruct(res)
   m14(:,1+i)=res(k3).RHO2;
   filename=strcat('Output/Figures/CSV/',model.filename,'_rho14.csv');
   writematrix(m14,filename,'Delimiter',';')
   system(['exec sed -i "s/\./,/g" ',filename]);
  end
  if ismember(15,plotfig)
   y4=real(model.fullEV(k3,1:numel(xPlot15)));
   m15(:,1+i)=y4;
   filename=strcat('Output/Figures/CSV/',model.filename,'_Lam15.csv');
   writematrix(m15,filename,'Delimiter',';')
  end
  if ismember(19,plotfig)
   y4=imag(model.fullEV(k3,1:numel(xPlot15)));
   m19(:,1+i)=y4;
   filename=strcat('Output/Figures/CSV/',model.filename,'_Lam19.csv');
   writematrix(m19,filename,'Delimiter',';')
  end
  if ismember(20,plotfig)
   m20(:,1+i)=res(k3).rhoDach;
   filename=strcat('Output/Figures/CSV/',model.filename,'_rho3D_20.csv');
   writematrix(m20,filename,'Delimiter',';')
   system(['exec sed -i "s/\./,/g" ',filename]);
  end
  if ismember(30,plotfig)
   m30(:,1+i)=res(k3).RHO2;
   filename=strcat('Output/Figures/CSV/',model.filename,'_rho30.csv');
   writematrix(m30,filename,'Delimiter',';')
   %system(['exec sed -i "s/\./,/g" ',filename]);
  end
  if ismember(32,plotfig)
   m32(:,1+i)=res(k3).RXB;
   filename=strcat('Output/Figures/CSV/',model.filename,'_RxB32.csv');
   writematrix(m32,filename,'Delimiter',';')
  end
  if ismember(34,plotfig)
   m34(:,1+i)=res(k3).DrhopDs;
   filename=strcat('Output/Figures/CSV/',model.filename,'_tRxB34.csv');
   writematrix(m34,filename,'Delimiter',';')
  end
 end
 
end
