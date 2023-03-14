%#!
%university:TU Wien 
clear 
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
 

   


   MyColours=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;0, 0, 1;0, 0.5, 0;1, 0, 0;0, 0.75, 0.75;0.75, 0, 0.75;0.75, 0.75, 0;0.25, 0.25, 0.25;1 1 0;1 0 1;0 1 1; 0 1 0;0 0 0 ];
     
    tmp=get(0,'defaultfigureposition');
 FesterPosX =uint16(tmp(1));
 FesterPosY =uint16(tmp(2));
 XBreite0=tmp(3);
 YHohe0=tmp(4);
Faktor=1;
 XBreite=XBreite0*Faktor;
 YHohe=YHohe0*Faktor;

 close all
   h314159=figure(314159);
  %set(gca,'Position',[0.10 0.08 0.90 0.92])%set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  %set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape', 'PaperType', 'A4')%,'Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   %set(h,'renderer','Painters')
   hold on

   wEVs=1:7;  
   %plot3(xMat2h(:,wEVs),yMat2h(:,wEVs),zMat2h(:,wEVs),'LineWidth',6,'Color','w')

   %wEVs=[1:3,6];
   myEVs=7;
   if exist('myEVs','var')
    if myEVs<7
     wEVs=myEVs;
     %wEVs=2
     for i=wEVs
      plot3(xMat2(:,i),yMat2(:,i),zMat2(:,i),'LineWidth',3,'Color',MyColours(i,:))
     end
    else
     wEVs=7;
     plot3(xMat2(:,wEVs),yMat2(:,wEVs),zMat2(:,wEVs),'LineWidth',3,'Color',MyColours(7,:),'LineStyle',':')
    end
    Exname=strcat(GFolder,'Figures/',name,'_3DIso',num2str(myEVs));
   else
    wEVs=1:7;
    %wEVs=2
    for i=wEVs
     plot3(xMat2(:,i),yMat2(:,i),zMat2(:,i),'LineWidth',3,'Color',MyColours(i,:))
    end
    wEVs=7;
    plot3(xMat2(:,wEVs),yMat2(:,wEVs),zMat2(:,wEVs),'LineWidth',3,'Color',MyColours(7,:),'LineStyle',':')
    Exname=strcat(GFolder,'Figures/',name,'_3DIso');
   end
   %plot(PT0([1,end],1),[0,0],'Color',[0 0 0],'LineWidth',1)
   xlabel('$\lambda$','Interpreter','latex')
   ylabel('$Re(\chi)$','Interpreter','latex')
   zlabel('$Im(\chi)$','Interpreter','latex')
   grid on
   %grid minor
   %set(gcf,'renderer','Painters')
   xlim([0 2.5])
   ylim([-1.5 1])
   zlim([-.5 .5])
   daspect([1 1 1])

   view([1,1,1])
   camroll(120)
   
   %exportgraphics(h,strcat(Exname,'.eps'),'BackgroundColor','none','ContentType','vector')
   print ('-dsvg', '-vector',strcat(Exname,'.svg'))
   %print('-dpng',strcat(Exname,'.png'))
   %print('-dpng','-r300',strcat(Exname,'300.png'))

  % fill(xMat2(1:424,1),yMat2(1:424,1),'r')
   %plot(xMat2(1:210,2),yMat2(1:210,2),'r')

