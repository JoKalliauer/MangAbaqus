%#!/bin/rm
%university:TU Wien 
 %filename='~/ownCloud/Post/MangAbaqus/Output/Figures/CSV/fein.csv';
 %filename='/home/jokalliau/ownCloud/Mang/Dokumente/own/DurchbiegungslinieB.csv';
 %filename='/home/jokalliau/mydataset.csv';
 %filename='/home/jokalliau/ownCloud/Mang/Dokumente/own/DurchbiegungslinieB_EnergieRho.csv';
 %filename='~/ownCloud/Post/MangAbaqus/Output/Figures/CSV/ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_rho14_JK.csv';
 %filename='/home/jokalliau/Abaqus/abaqus.rpt';
 %filename='/home/jokalliau/Abaqus/Durchbeigung/both16.csv';
 %filename='/home/jokalliau/DurchbiegungMathematicaNSVgl.csv';
 filename='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.02-KNL2_rho3D_20.csv';
 
 %opts = detectImportOptions(filename,'FileType','fixedwidth')
 opts = delimitedTextImportOptions(...
'VariableNames',{'lambda' 'rho1st' 'rho2nd' 'rho3nd'},'DataLine',[1 Inf],...
'NumVariables',4,'Delimiter',{';','\t'}); %'VariableWidths',[6 12 12 12 12]
system(['exec sed -i "s/,/\./g" ',filename])
 T0 = readtable(filename,opts);
 %PT0=str2double(table2array(T0));
 PT0=Commastr2doubleJK(table2array(T0));
 %PT0=table2array(T0);
%  outputfile='~/Abaqus/Matlab.csv';
%  writematrix(PT0,outputfile,'Delimiter',';')
%  system(['exec sed -i "s/\./,/g" ',outputfile]);
 
 close all
  figure(14);
  set(gca,'Position',[0.10 0.08 0.90 0.92])%set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  %set(gca,'xdir', 'reverse')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape', 'PaperType', 'A4')%,'Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on

   grid on
   grid minor

  xlabel('Lamdba $\lambda$','Interpreter','latex');
  %ylabel('$w_m$','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  %title('deltalambda=0.005')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %daspect([1 1 1])
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(PT0(:,1),PT0(:,2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',[0, 0.4470, 0.7410]
  %legend('Energieverhaelntis Abaqus','Energieverhaeltnis Th. II. Ord.','Kruemmungsradius Abaqus','Location','southeast')
  %yticks(0:.05:1)
  %bbb.XLim = [0 2.2];

   print('-dsvg',strcat(filename,'.svg'))
   print('-dpng',strcat(filename,'.png'))
   print('-dpdf',strcat(filename,'.pdf'),'-fillpage')

