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
GFolder='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/'; ylabelJK='radius of the first Frenet-curvature $\rho_1$'; %%#ok<NASGU>
%diagramname='TL_arch3D-B32-20-loadfac-1-eps0.005-KNL2_rho30.csv';xlabelJK='line load $p$ [$\textrm{N}/\textrm{m}$]'; 
 %buck=2.77E6; cfig('yLim')=[0 1]; cfig('xLim')=[0 6E6]; cfig('Cxticks')=[0:1E6:2E6 buck 3E6:1E6:6E6]; cfig('xticklabels')={0 '1E6' '2E6' '$p_S$~~' '~~3E6' '4E6' '5E6' '6E6'}; cfig('xline')=buck;
%diagramname='TL_arch3D-B32-20-loadfac-1-eps0.005-I_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$'; xlabelJK='line load $p$ [$\textrm{N}/\textrm{m}$]';
 % buck=2.77E6;  cfig('yLim')=[0 1]; cfig('xLim')=[0 6E6]; cfig('xline')=buck; cfig('Cxticks')=[0:1E6:2E6 buck 3E6:1E6:6E6]; cfig('xticklabels')={0 '1E6' '2E6' '$p_S$~~' '~~3E6' '4E6' '5E6' '6E6'};
%diagramname='pureBendingBeamJK-B32OSH-20-len-5-loadfac-1-eps0.02-KNL2_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M$ [N\,m]'; buck=286E3;
 %cfig('xLim')=[0 2E6]; cfig('Cxticks')=[0 buck 0.5E6:0.5E6:2E6 3E6:1E6:6E6]; cfig('xticklabels')={0 '$M_S$' '5E5' '10E5' '15E5' '20E5' '3E6' '4E6' '5E6' '6E6'};cfig('xline')=buck;
%diagramname='pureBendingBeamJK-B32OSH-20-len-5-loadfac-1-eps0.02-I_rho30.csv'; ylabelJK='radius of the first Frenet-curvature $\rho_1$';xlabelJK='bending moment $M$ [N\,m]'; buck=286E3;
 %cfig('xLim')=[0 2E6]; cfig('Cxticks')=[0 buck 0.5E6:0.5E6:2E6 3E6:1E6:6E6]; cfig('xticklabels')={0 '$M_S$' '5E5' '10E5' '15E5' '20E5' '3E6' '4E6' '5E6' '6E6'};cfig('xline')=buck;
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_rho30.csv';xlabelJK='normal force $N$ [N]'; buck=Inf;
 %cfig('xLim')=[0 1.5E6]; cfig('Cxticks')=[0 0.5E6:0.5E6:2E6 3E6:1E6:6E6]; cfig('xticklabels')={0 '5E5' '10E5' '15E5' '20E5' '3E6' '4E6' '5E6' '6E6'};
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-KNL2_RxB32.csv'; xlabelJK='normal force $N$ [N]'; buck=Inf; cfig('yLim')=[-1 1]; ylabelJK='$\mathbf{r}_1\cdot\mathbf{e}_3$';
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-I_rho30_edit.csv';xlabelJK='normal force $N$ [N]'; buck=Inf;
%diagramname='ecc-B32OS-20-len-5-ecc-0.040447-loadfac-1-eps0.01-I_RxB32_edit.csv';xlabelJK='normal force $N$ [N]'; buck=Inf; cfig('yLim')=[-1 1]; ylabelJK='$\mathbf{r}_1\cdot\mathbf{e}_3$';

%buck=2.77E6; cfig('xLim')=[1E6 6E6]; cfig('Cxticks')=[0:1E6:2E6 buck 3E6:1E6:6E6]; cfig('xticklabels')={0 '1E6' '2E6' '$p_S$~~' '~~3E6' '4E6' '5E6' '6E6'}; cfig('xline')=buck; 
xlabelJK='line load $p$ [$\textrm{N}/\textrm{m}$]';
cfig('xLim')=[.2 .7]; xlabelJK='Amplification factor $\lambda$ [-]';

% name='TL15';
% diagramnamei{1}='TL_arch3D-B32H-20-loadfac-1-eps0005-KNL2_tRxB34_edit.csv';
% diagramnamei{2}='TL_arch3D-B32H-20-loadfac-1-eps001-KNL2_tRxB34_edit.csv';
% diagramnamei{3}='TL_arch3D-B32H-20-loadfac-1-eps002-KNL2_tRxB34_edit.csv';
% ylabelJK='$\kappa_2\,(\mathbf{r}_1\cdot\mathbf{e}_3)$';
% % %cfig('yLim')=[-30 20];
% cfig('yLim')=[-10 15];

% name='TLeb';
% diagramnamei{1}='TL_arch3D-B32H-20-loadfac-1-eps0005-KNL2_RxB32_edit.csv';
% diagramnamei{2}='TL_arch3D-B32H-20-loadfac-1-eps001-KNL2_RxB32_edit.csv';
% diagramnamei{3}='TL_arch3D-B32H-20-loadfac-1-eps002-KNL2_RxB32_edit.csv';
% ylabelJK='$\mathbf{r}_1\cdot\mathbf{e}_3$';

name='TLtau';
diagramnamei{1}='TL_arch3D-B32H-20-loadfac-1-eps0005-KNL2_tau07_edit.csv';
diagramnamei{2}='TL_arch3D-B32H-20-loadfac-1-eps001-KNL2_tau07_edit.csv';
diagramnamei{3}='TL_arch3D-B32H-20-loadfac-1-eps002-KNL2_tau07_edit.csv';
ylabelJK='$\kappa_2$';
% %cfig('yLim')=[-30 20];
% %cfig('yLim')=[-15 10];



% opts = detectImportOptions(filename,'FileType','fixedwidth');
 opts = delimitedTextImportOptions(...
'VariableNames',{'lambda' 'rho1st' 'rho2nd' 'rho3nd'},'DataLine',[1 Inf],...
'NumVariables',4,'Delimiter',{';','\t'}); %'VariableWidths',[6 12 12 12 12]
%system(['exec sed -i "s/,/\./g" ',filename])


 PT01=Commastr2doubleJK(table2array(readtable(strcat(GFolder,diagramnamei{1}),opts)));
 PT02=Commastr2doubleJK(table2array(readtable(strcat(GFolder,diagramnamei{2}),opts)));
 PT03=Commastr2doubleJK(table2array(readtable(strcat(GFolder,diagramnamei{3}),opts)));
 
 xVal=unique([PT01(:,1);PT02(:,1);PT03(:,1)]);
 len=size(xVal,1);
 yMat=NaN(len,3);
 
  Pos1=NaN(size(PT01,1),1);
 for i=1:size(PT01,1)
  Pos1(i)=find(xVal==PT01(i));
 end
 yMat(Pos1,1)=PT01(:,2);
 
 Pos2=NaN(size(PT02,1),1);
 for i=1:size(PT02,1)
  Pos2(i)=find(xVal==PT02(i));
 end
 yMat(Pos2,2)=PT02(:,2);
 
  Pos3=NaN(size(PT03,1),1);
 for i=1:size(PT03,1)
  Pos3(i)=find(xVal==PT03(i));
 end
 yMat(Pos3,3)=PT03(:,2);

 if strcmp(name,'TL15')
  yMat=-yMat;
 end
 
 ltp=8330000;
 if strcmp(name,'TL15')|| strcmp(name,'TLtau')
  PT0=[xVal,yMat];
 else
  PT0=[xVal/ltp,yMat];
 end

diagramname='TL';

lenJK=size(PT0,1);
 
% close all
  figure(14);
  set(gca,'Position',[0.10 0.08 0.90 0.92])%set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  %set(gca,'xdir', 'reverse')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape', 'PaperType', 'A4')%,'Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on

   grid on
   grid minor
   
   %xlabelJK='Lamdba $\lambda$';
  xlabel(xlabelJK,'Interpreter','latex');
  %ylabel('$w_m$','Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %title('deltalambda=0.005')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %daspect([1 1 1])
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(PT0(:,1),PT0(:,2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',[0, 0.4470, 0.7410]
  plot(0,PT0(5,2),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  text(0,PT0(5,2),'$\rho_0$','interpreter','latex')
  %legend('Energieverhaelntis Abaqus','Energieverhaeltnis Th. II. Ord.','Kruemmungsradius Abaqus','Location','southeast')
  %yticks(0:.05:1)
  %bbb.XLim = [0 2.2];

   
   
   
   diagramname(diagramname == ' ' | diagramname == '/'| diagramname == '.'| diagramname == '_') = [];
   
   cfig('order')=NaN;
   %cfig('lineColor')={[1 0 0],[1 0 0]};
   cfig('LineStyle')={'-','--',':'};
   cfig('MyMarker')='none';
   cfig('MyMarkerMulti')={'none' 'none' 'none'};
   %cfig('xLim')=[0 1.5E6];
   cfig('wofuer')='TeX';
   %cfig('Cxticks')=[0 buck 1E6:1E6:6E6];
   %cfig('xticklabels')={0 '$M_S$' '1E6' '2E6' '3E6' '4E6' '5E6' '6E6'};
   cfig('LineWidth')=3.5;
   cfig('XAxisLocation')='bottom';
   cfig('lineColor')='JKextended';
   cfig('closing')=NaN;

   plotitJK(PT0(:,1),yMat,GFolder,xlabelJK,ylabelJK,diagramname,cfig)
   grid on
   grid minor
   plot(PT0([1,end],1),[0,0],'Color',[0 0 0],'LineWidth',1.5)
   lh=legend('$\Delta\lambda=0.005$','$\Delta\lambda=0.01$','$\Delta\lambda=0.02$','Interpreter','latex','Location','northeast');
%    if strcmp(name,'TL15')
%     lh.Location='southeast';
%    end
   print ('-depsc',strcat(GFolder,'Figures/',name,'.eps'))
   
   diagramname %#ok<NOPTS>
   

