%#!
%university:TU Wien 
cfig = containers.Map;
GFolder='/home/jkalliau/ownCloud/Post/MangAbaqus/Output/Figures/CSV/'; %%#ok<NASGU>
diagramname='EccEW.csv';
ylabelJK='eigenvalue $\chi$'; 
close all
 
cfig('xLim')=[0 2.5]; xlabelJK='Amplification factor $\lambda$ [-]';

name='EW22';
 opts = delimitedTextImportOptions(...
'VariableNames',{'lambda' 'rho1st' 'rho2nd' 'rho3nd'},'DataLine',[2 Inf],...
'NumVariables',4,'Delimiter',{';','\t'}); %'VariableWidths',[6 12 12 12 12]
%system(['exec sed -i "s/,/\./g" ',filename])


 PT01=Commastr2doubleJK(table2array(readtable(strcat(GFolder,diagramname),opts)));
 
 xVal=unique([PT01(:,1)]);
 xValCor=xVal./2.2;
 len=size(xVal,1);
 yMat=NaN(len,6);
 
  Pos1=NaN(size(PT01,1),1);
 for i=1:size(PT01,1)
  Pos1(i)=find(xVal==PT01(i));
 end
 yMat(Pos1,1:8)=PT01(:,[2,3,6,4,4,7,8,9]);
 zMat=0*yMat;
 %xMat=zMat;
 zMat(Pos1,4)=PT01(:,5);
 zMat(Pos1,5)=-PT01(:,5);

 xMat=[xValCor,xValCor,xValCor,xValCor,xValCor,xValCor,xValCor,xValCor];

 xy2=[xMat(:,2),yMat(:,2)];
 xy2d = sortrows(xy2,2,'descend','MissingPlacement','last');
 xMat2=xMat;
 yMat2=yMat;
 zMat2=zMat;

 xMat2(:,2)=xy2d(:,1);
 yMat2(:,2)=xy2d(:,2);
 zMat2(1:424,2)=zeros(424,1);
 zMat2h=zMat2-.01;
 yMat2h=yMat2-.001;
 xMat2h=xMat2-.001;



  PT0=[xValCor,yMat];

diagramname='old';

lenJK=size(PT0,1);
 
% close all
  figure(14);
  set(gca,'Position',[0.10 0.08 0.90 0.92])%set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape', 'PaperType', 'A4')%,'Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on

   grid on
   grid minor
   
  xlabel(xlabelJK,'Interpreter','latex');
  bbb = gca();
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xMat2,yMat2,'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',[0, 0.4470, 0.7410]

   
   
   
   diagramname(diagramname == ' ' | diagramname == '/'| diagramname == '.'| diagramname == '_') = [];
   
   cfig('order')=NaN;
   %cfig('LineStyle')={'-','--',':'};
   cfig('MyMarker')='none';
   %cfig('MyMarkerMulti')={'none' 'none' 'none'};
   cfig('wofuer')='TeX';
   cfig('LineWidth')=3.5;
   cfig('XAxisLocation')='bottom';
   cfig('lineColor')='JKextended';
   cfig('closing')=NaN;

   %plotitJK(PT0(:,1),yMat,GFolder,xlabelJK,ylabelJK,diagramname,cfig)
   
   grid on
   grid minor
   plot(PT0([1,end],1),[0,0],'Color',[0 0 0],'LineWidth',1.5)
   %lh=legend('$\Delta\lambda=0.005$','$\Delta\lambda=0.01$','$\Delta\lambda=0.02$','Interpreter','latex','Location','northeast');
   print ('-depsc',strcat(GFolder,'Figures/',name,'.eps'))
   


   MyColours=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;0, 0, 1;0, 0.5, 0;1, 0, 0;0, 0.75, 0.75;0.75, 0, 0.75;0.75, 0.75, 0;0.25, 0.25, 0.25;1 1 0;1 0 1;0 1 1; 0 1 0;0 0 0 ];
     
 close all
   h=figure(314159);
   %set(h,'renderer','Painters')
   hold on

   wEVs=1:7;
   plot3(xMat2h(:,wEVs),yMat2h(:,wEVs),zMat2h(:,wEVs),'LineWidth',6,'Color','w')

   wEVs=1:6;
   for i=wEVs
    plot3(xMat2(:,i),yMat2(:,i),zMat2(:,i),'LineWidth',3,'Color',MyColours(i,:))
   end
   wEVs=7;
   plot3(xMat2(:,wEVs),yMat2(:,wEVs),zMat2(:,wEVs),'LineWidth',3,'Color',MyColours(7,:),'LineStyle',':')
   plot(PT0([1,end],1),[0,0],'Color',[0 0 0],'LineWidth',1)
   xlabel('$\lambda$','Interpreter','latex')
   ylabel('$Re(\chi)$','Interpreter','latex')
   zlabel('$Im(\chi)$','Interpreter','latex')
   grid on
   %grid minor
   xlim([0 2.5])
   daspect([1 1 1])

   view([1,1,1])
   camroll(120)
   Exname=strcat(GFolder,'Figures/',name,'_3DIso.eps');
   exportgraphics(h,Exname,'BackgroundColor','none','ContentType','vector')
   print('-dpng','Iso.png')

   view([0,0,1])
   print ('-depsc',strcat(GFolder,'Figures/',name,'_3Dxy.eps'))
   print('-dpng','xy.png')

   view([0,1,0])
   camroll(180)
   print -depsc -tiff -r300 -vector 3Dxz.eps
   print('-dpng','xz.png')

   view([1,0,0])
   camroll(90)
exportgraphics(h,'3Dxyz.eps','BackgroundColor','none','ContentType','vector')
   print('-dpng','yz.png')

