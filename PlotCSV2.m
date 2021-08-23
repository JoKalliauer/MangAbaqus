%#!
%university:TU Wien 
%filename='~/ownCloud/Post/MangAbaqus/Output/Figures/CSV/fein.csv';
%filename='/home/jokalliau/ownCloud/Mang/Dokumente/own/DurchbiegungslinieB.csv';
%filename='/home/jokalliau/mydataset.csv';
%filename='/home/jokalliau/ownCloud/Mang/Dokumente/own/DurchbiegungslinieB_EnergieRho.csv';
%filename='~/ownCloud/Post/MangAbaqus/Output/Figures/CSV/ecc-B32OS-100-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_rho14_JK.csv';
%filename='/home/jokalliau/Abaqus/abaqus.rpt';
%filename='/home/jokalliau/Abaqus/Durchbeigung/both16.csv';
%filename='/home/jokalliau/DurchbiegungMathematicaNSVgl.csv';
%filename='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/ecc-B32OS-2-len-5-ecc-0.16467-loadfac-1-eps0.02-KNL2_rho3D_20.csv';
%filename='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/TL_arch3D-B32-20-loadfac-1-eps0.005-KNL2_rho30.csv';
cfig = containers.Map;
GFolder='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/'; ylabelJK='radius of the first Frenet-curvature $\rho_1$'; Faktor=1.54; %#ok<*NASGU> 
%diagramname='TL_arch3D-B32-20-f1-eps0.005-u1-KNL2_rho30.csv';xlabelJK='line load $p$ [$\textrm{N}/\textrm{m}$]'; name='Fig2a'
 %buck=2.77E6; cfig('yLim')=[0 1]; cfig('xLim')=[0 3E6]; %cfig('Cxticks')=[0:1E6:2E6 buck 3E6:1E6:6E6]; cfig('xticklabels')={0 '1E6' '2E6' '$p_S$~~' '~~3E6' '4E6' '5E6' '6E6'};
%diagramname='TL_arch3D-B32-20-f1-eps0.005-u1-I_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$'; xlabelJK='line load $p$ [$\textrm{N}/\textrm{m}$]'; name='Fig2b';
 %buck=2.77E6;  cfig('yLim')=[0 1]; cfig('xLim')=[0 3E6]; %cfig('Cxticks')=[0:1E6:2E6 buck 3E6:1E6:6E6]; cfig('xticklabels')={0 '1E6' '2E6' '$p_S$~~' '~~3E6' '4E6' '5E6' '6E6'};
%alt %diagramname='pureBendingBeamJK-B32OSH-20-len-5-loadfac-1-eps0.02-KNL2_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M_y$ [N\,m]'; buck=286E3; name='Fig4a';
%diagramname='BB5-B32OSH20-l5-f1-eps0.02-u1-KNL2_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M_y$ [N\,m]'; buck=286E3; name='Fig4a';
cfig('xLim')=[0 3E5];% cfig('Cxticks')=[0 buck 0.5E6:0.5E6:2E6 3E6:1E6:6E6]; cfig('xticklabels')={0 '$(M_y)_S$~' '~~5E5' '10E5' '15E5' '20E5' '3E6' '4E6' '5E6' '6E6'};
%alt %diagramname='pureBendingBeamJK-B32OSH-20-len-5-loadfac-1-eps0.02-I_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M_y$ [N\,m]'; buck=286E3;name='Fig4b';
diagramname='BB5-B32OSH20-l5-f1-eps0.005-u1-I_rho30_edit.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M_y$ [N\,m]'; buck=286E3;name='Fig4b';

%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_rho30.csv'; buck=Inf;name='Fig7a';
 %cfig('xLim')=[0 1.5E6]; cfig('Cxticks')=[0 0.5E6:0.5E6:2E6 3E6:1E6:6E6]; cfig('xticklabels')={0 '500' '1000' '1500' '20E2' '3E3' '4E3' '5E3' '6E3'};xlabelJK='applied force $P$ [kN]'; %Faktor=1/1000;
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_RxB32.csv'; buck=Inf; cfig('yLim')=[-1 1]; ylabelJK='$\mathbf{r}_1\cdot\mathbf{e}_3$';name='Fig8a';
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-I_rho30_edit.csv'; buck=Inf;name='Fig7b';
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-I_RxB32_edit.csv'; buck=Inf; cfig('yLim')=[-1 1]; ylabelJK='$\mathbf{r}_1\cdot\mathbf{e}_3$';name='Fig8b';



 filename=strcat(GFolder,diagramname);
 
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

lenJK=size(PT0,1);
 
 close all
  figure(14);
  set(gca,'Position',[0.10 0.08 0.90 0.92])%set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  %set(gca,'xdir', 'reverse')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape', 'PaperType', 'A4')%,'Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on

   %grid on
   %grid minor
   
   %xlabelJK='Lamdba $\lambda$';
  xlabel(xlabelJK,'Interpreter','latex');
  %ylabel('$w_m$','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  %title('deltalambda=0.005')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %daspect([1 1 1])
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(PT0(:,1),PT0(:,2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',[1 0 0]);%,'Color',[0, 0.4470, 0.7410]
  plot(0,PT0(5,2),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  text(0,PT0(5,2),'$\rho_0$','interpreter','latex')
  %legend('Energieverhaelntis Abaqus','Energieverhaeltnis Th. II. Ord.','Kruemmungsradius Abaqus','Location','southeast')
  %yticks(0:.05:1)
  %bbb.XLim = [0 2.2];

   print('-dsvg',strcat(filename,'.svg'))
   print('-dpng',strcat(filename,'.png'))
   print('-dpdf',strcat(filename,'.pdf'),'-fillpage')
   
   
   
   diagramname(diagramname == ' ' | diagramname == '/'| diagramname == '.'| diagramname == '_') = [];
   
   cfig('order')=NaN;
   cfig('lineColor')={[1 0 0],[1 0 0]};
   cfig('LineStyle')={'-','-'};
   cfig('MyMarker')='none';
   cfig('MyMarkerMulti')={'none' 'none'};
   %cfig('xLim')=[0 1.5E6];
   cfig('wofuer')='TeX';
   %cfig('Cxticks')=[0 buck 1E6:1E6:6E6];
   %cfig('xticklabels')={0 '$M_S$' '1E6' '2E6' '3E6' '4E6' '5E6' '6E6'};
   cfig('LineWidth')=5;
   cfig('XAxisLocation')='bottom';
   cfig('grid')='off';
   cfig('Faktor')=Faktor;
   PT0a=NaN(lenJK,1);
   PT0b=NaN(lenJK,1);
   PT0a(PT0(:,1)<buck)=PT0(PT0(:,1)<buck,2);
   PT0b(PT0(:,1)>buck)=PT0(PT0(:,1)>buck,2);
   plotitJK(PT0(:,1),PT0(:,2),GFolder,xlabelJK,ylabelJK,diagramname,cfig)
   if strcmp(diagramname,'TLarch3D-B32-20-loadfac-1-eps0005-KNL2rho30csv') || strcmp(diagramname,'TLarch3D-B32-20-loadfac-1-eps0005-Irho30csv') || strcmp(diagramname,'pureBendingBeamJK-B32OSH-20-len-5-loadfac-1-eps002-KNL2rho30csv')
    plot([buck buck],[0,PT0a(find(PT0(:,1)<buck,1,'last'))],'Color',[0 0 0],'LineWidth',2)
   end
   if strcmp(diagramname,'ecc-B32OS-20-len-5-ecc-0040447-loadfac-1-eps001-KNL2RxB32csv') ||strcmp(diagramname,'ecc-B32OS-20-len-5-ecc-0040447-loadfac-1-eps001-IRxB32editcsv')
    plot(PT0([1,end],1),[0,0],'Color',[0 0 0],'LineWidth',2)
   end
   if strcmp(diagramname,'ecc') || strcmp(diagramname,'ecc-B32OS-20-len-5-ecc-0040447-loadfac-1-eps001-KNL2rho30csv')
    plot(0,PT0(5,2),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
    text(0,PT0(5,2)+.03,'\,$ (\rho_1)_0$','interpreter','latex','FontSize',25)
   end
   if strcmp(diagramname,'ecc-B32OS-20-len-5-ecc-0040447-loadfac-1-eps001-Irho30editcsv')
    plot(0,PT0(5,2),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
    text(0,PT0(5,2)-.03,'\,$ (\rho_1)_0$','interpreter','latex','FontSize',25)
   end
   print ('-depsc',strcat(GFolder,'Figures/',name,'.eps'))
   
   diagramname %#ok<NOPTS>
   

