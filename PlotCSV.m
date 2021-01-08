%#!/bin/rm
%university:TU Wien 
 %filename='~/ownCloud/Post/MangAbaqus/Output/Figures/CSV/fein.csv';
 filename='/home/jokalliau/ownCloud/Mang/Dokumente/own/DurchbiegungslinieB.csv';
 %opts = detectImportOptions(filename)
 opts = delimitedTextImportOptions(...
'VariableNames',{'lambda' 'rho1st' 'rho2nd' 'rho3nd'},'DataLine',[1 Inf],...
'NumVariables',4,'Delimiter',';'); %'VariableWidths',[6 12 12 12 12]
 T0 = readtable(filename,opts);
 PT0=Commastr2doubleJK(table2array(T0));
 
 close all
  figure(14);
  %set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on

   grid on
   grid minor

  xlabel('lambda');
  ylabel('$\rho$','Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0.0,1];
  title('deltalambda=0.005')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(PT0(:,1),PT0(:,2:end),'LineStyle','-','Marker','none','LineWidth',1.5);
  %yticks(0:.05:1)
  bbb.XLim = [0 inf];

   print('-dsvg',strcat(filename,'_rho14.svg'))
   print('-dpng',strcat(filename,'_rho14.png'))
   print('-dpdf',strcat(filename,'_rho14.pdf'),'-fillpage')

