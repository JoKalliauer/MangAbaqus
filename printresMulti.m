function printresMulti(res,model,plotfig,~,~,resEWs,~)


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
  

 k3=resEWs;
 %lambda = res(1).lambda;
 m14=NaN(numel(res(1).lambda),1+numel(k3));
 m14(:,1)=res(1).lambda;
 m20=m14;
 for k3=resEWs
  i=find(resEWs==k3);
  if ismember(14,plotfig)
   m14(:,1+i)=res(k3).RHO2;
   filename=strcat('Output/Figures/CSV/',model.filename,'_rho14.csv');
   writematrix(m14,filename,'Delimiter',';')
   system(['exec sed -i "s/\./,/g" ',filename]);
  end
  if ismember(20,plotfig)
   m20(:,1+i)=res(k3).rhoDach;
   filename=strcat('Output/Figures/CSV/',model.filename,'_rho3D_20.csv');
   writematrix(m20,filename,'Delimiter',';')
   system(['exec sed -i "s/\./,/g" ',filename]);
  end
 end
 
end
