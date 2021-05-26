 function plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main)
 % 2...rho
 % 1...chi
 % 3..._velocity
 % 4_tanacceleration
 % 5_noracceleration
 % 6_totacceleration
 % 7_torque
 % 8_SinCos
 % 8.1 Sin/Cos
 % 11_ Eigenwert ylim([0,1])
 % % 11.1 ylim detail
 % 12_real(Eigenwert) kein ylim:-
 % 13_abs(eigenwert)
 % 14_rho2
 % 15_real EW..fulllambda
 % 16 d^2 EW / (dlambda)^2
 %  16.1 only negativ values
 % 17 OC6
 % 18 OC7
 % 19 imag(EW)
 % 20 rhorho3D,rhoDach
 % 21 singamma
 % 211 debugging
 % 22 OC5
 % 23 Hypo
 % 24 =>7
 % 25 CosGamma
 % 26 r*b
 % 27 r*b (3D)
 %28 Rconst
 %29 EBENE
 %30 rho  -- p
 %31 1-rho -- p
 %32 r*b -- p
 %33 arccos(r-r(l))
 %34 DrhopDs
 %35 chi -- p
 %36 eigenwert detail
 %37 Eigenwert nonlinear
 %38 cosPhiMang
 %39 ZaelerB
 %40 NennerB1
 %41 NennerB2
 %42 NormR
 
 
 %900 detKt
 %901 sqrt[N]{detKt}
 %902 \xiJK
 %906 dl/dxi
 %908 v [m]
 %909 d\lambda/dw
 %914 det/max(det)
 
 if numel(resEWs)>0
  lambda = res(min(resEWs)).lambda;
 else
  lambda=model.lambda;
  %markJK=MyMarker(mod(main.colorshift-1,2)+1);
  %markJK='none';
  %colJK=MyColours{mod(main.colorshift-1,19)+1};
 end
 lengthlam=length(lambda);
 lamlast=find(~isnan(lambda),1,'last');

 export=true;
 MyLines={'-','--','-.',':'};
   cfig = containers.Map;
   cfig('ylog')=false;
   cfig('yMin')=0;
   %cfig('yMax')=1;
   cfig('xMin')=0;
   cfig('line')='-';
   cfig('wofuer')='TeX';
   cfig('lineColor')=[0, 0.4470, 0.7410];
   cfig('MyMarker')='none';
   cfig('LineWidth')=3;
   cfig('FontSize')=17.5;
   cfig('order')=NaN;
   cfig('closing')=NaN;
   cfig('LineStyle')=MyLines{mod(main.colorshift,4)+1};
%    cfig('lineColor')=colJK;
 if export && isfield(model,'load')
  lengthlam=min(lengthlam,numel(model.load)+1);
  if strcmp(model.filename(1),'p') %%#ok<UNRCH>
   xlabelload='bending moment $M$'; % [kN\,m]
   xValload=model.load(1:lengthlam-1);%/1000;
   xValload0=model.load0;%/1000;
   cfig('xMax')=5000;
  elseif strcmp(model.filename(1),'T')
   xlabelload='line load $p$ [kN/cm]';
   xValload=model.load(1:lengthlam-1);%/100000;
   xValload0=model.load0;%/100000;
   %cfig('xMax')=70;
   %cfig('Cxticks')=0:10:70;
  elseif strcmp(model.filename(1),'e')
   xlabelload='normal force $N$ [N]';
   xValload=model.load(1:lengthlam-1);%/1000;
   xValload0=model.load0;%/1000;
   %cfig('xMax')=1600;
   %cfig('Cxticks')=0:200:1600;
  else
   xlabelload='load';
   xValload=model.load(1:lengthlam-1);
   xValload0=model.load0;
  end
%   xValload0=[0;xValload(1:end)];
  NrLoad0=numel(xValload0);
  xValfullload0=model.fullload0;
 end
 if sum(strcmp(fieldnames(model), 'xlabelloadname')) ~= 0
   xlabelload=model.xlabelloadname;
 end
 
 
 
  %% Figuredef
  FontSize0=(get(0,'defaultaxesfontsize')+get(0,'defaulttextfontsize'))/2;
  screenX=1920;
  %screenY=1080;
 tmp=get(0,'defaultfigureposition');
 %FesterPosX =uint16(tmp(1));
 FesterPosY =uint16(tmp(2));
 XBreite0=tmp(3);
 YHohe0=tmp(4);

FontName = 'times';
Faktor=1;
XBreite=(XBreite0)*Faktor; %416
YHohe=(YHohe0)*Faktor; %312
FontSize=FontSize0*Faktor;

%gcfPostition=[FesterPosX   FesterPosY   XBreite   YHohe]; 
FesterPosXNR=uint16(linspace(0,screenX-XBreite,numel(plotfig)));

%% 

 
 if ~exist('model','class')
  if isfield(model,'filename')
   modelfilename=model.filename;
  else
   modelfilename=model.TestCase;
  end
 else
  modelfilename='unknown';
 end
 if ~exist('plotfig','var')
  plotfig=[1,2];
 end
 %MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
 
 
  if ismember(0,plotfig)
  figure(2147483646);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==0)   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='lambda';
  ylabel(ylabelJK,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  Points=numel(xValload);
  plot(xValload,lambda(1:Points),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',MyColours{1});
 end

 
 
 %maxlam = max(lambda);
 %maxlam = 1.0;
 if ~exist('main','var')
  main.savefigures=false;
  main.colorshift=0;
 end
 if main.savefigures==true
  if ~exist('Output/Figures/', 'dir')
   if isunix
    mkdir('Output/Figures/SVG');
    mkdir('Output/Figures/PNG');
    mkdir('Output/Figures/PDF');
   elseif ispc
    warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
    return
   elseif ismac
    warning('MyProgram:OS','Mac not tested try unix')
   end
  end
 end
 model.savefigures=false;%might get overwritten for k3==resEWs(end)
 %set(0, 'DefaultFigureWindowState', 'minimized');
 cfig('Faktor')=1;
 for k3=resEWs
  %k3 %#ok<NOPRT>
  if k3==resEWs(end)
   model.savefigures=main.savefigures;
  end
  Nr=k3+main.colorshift-1;
  colJK=MyColours{mod(Nr,19)+1};
  cfig('lineColor')=colJK;
%   lineStyleJK=MyLines{mod(Nr,4)+1};
  %lambda1 = res.lambdaorg;
  lambda = res(k3).lambda;
  if numel(lambda)<=41
   markJK=MyMarker(mod(Nr,2)+1);
  else
   markJK='none';
  end
% markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
 %V = res.V;
 %A = res.A;
 S = res(k3).S;
 A0 = res(k3).A0;%...total accerlation
 At = res(k3).At;
 An = res(k3).An;
 RHO = res(k3).RHO;
 %RHO2 = res.RHO2;
 %TAU = res(k3).TAU;
 cosmu = res(k3).cosmu;
 sinpsi = res(k3).sinpsi;
 %       OC1 = OC1(Orth);
 %OCeig = res.OCeig;
 LAM = res(k3).LAM;
 %fulllambda=model.fulllambda
 %fullEV=real(model.fullEV(k3,:));
 %R = res.R;
 %POS = res.POS;
 
 %OC0 = res.OC0;
 %OC1 = res.OC1;
 %OC2 = res.OC2;
 %OC3 = res.OC3;
 
 %coplanar = res.coplanar;
 
 %colo = [rand(), rand(), rand()];
 %colo = [0 0 0];
 
 %kl = res.stability_limit(1); 
 %stability_limit = res.stability_limit(2);
 
 

 %lambda = lambda/maxlam;
 

 

 
 lambda = reshape(lambda,numel(lambda),1);
 LAM = reshape(LAM,length(LAM),1);
 %LAM = LAM/maxlam;
 if ismember(1,plotfig)
  figure(1);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  %        plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none');%,'Color',colo);
  %        plot(lambda(2:end), LAM(2:end),'LineStyle','--','Marker','none');%,'Color',colo);
  aaa = gca();
  if size(aaa.Children,1)<2
   plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
   plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
  end
  plot(lambda(2:end), lambda(2:end)+real(LAM(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);% ,'Color',col %,'Color',colo);
  xlabel('Lambda');
  %        ylabel('Lambda + chi | or | chi');
  ylabel('EW');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %        legend('chi','eigenvalue','lambda','zero')
  title(modelfilename)
  %xticks(0:10)
  %yticks(-40:50)
  %ylim([0 min(20,lambda(end)+real(LAM(end)))])
  %daspect([1 1 1])
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_chi.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_chi.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_chi.pdf'),'-dpdf')
  end
 end

 if ismember(2,plotfig)
  figure(2);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  plot(lambda(2:end),RHO(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  
  %        plot(lambda(2:end),RHO2(2:end),'LineStyle','none','Marker','.','Color',colo);
  %        plot(lambda,RHO3,'LineStyle','none','Marker','o','Color',colo);
  %        plot(lambda,RHO4,'LineStyle','--','Marker','x','Color',colo);
  xlabel('lambda');
  ylabel('rho');
  %title('rho');
  bbb = gca();
  %bbb.YLim = [0.0,1.1];
  bbb.XLim = [0 inf];
  %col = bbb.Children(1).Color;
  
  %        plot(lambda(kl),RHO(kl),'LineStyle','none','Marker','x','LineWidth',1.5,'Color',col);
  grid on
  title(modelfilename,'Interpreter','none')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho.pdf'),'-dpdf')
  end
 end
 if ismember(-2,plotfig)
  disp([lambda(2:end),RHO(2:end)])
 end
 
 if ismember(3,plotfig)
  figure(3)
  %set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  plot(lambda(2:end),S(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('velocity');
  title('velocity');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %nopeaks = (S(S<10*mean(S(:))));
  %bbb = gca();
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_velocity.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_velocity.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_velocity.pdf'),'-fillpage')
  end
 end
 
 if ismember(3.3,plotfig)
  figure(33)
  %set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==3.3)   FesterPosY   XBreite   YHohe]); 
  hold on
  plot(lambda(2:end),res(k3).SpeedDach(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('velocity $\dot{\hat{s}}$','Interpreter','latex');
  title('velocity');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %nopeaks = (S(S<10*mean(S(:))));
  %bbb = gca();
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_velocity3D.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_velocity3D.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_velocity3D.pdf'),'-fillpage')
  end
 end
 
 if ismember(4,plotfig)
  figure(4)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  At = real(At);
  plot(lambda(2:end),At(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('tangential acceleration');
  title('tangential acceleration');
  grid on
  %nopeaks = (At(abs(At)<10*mean(abs(At(:)))));
  %bbb = gca();
  %bbb.YLim = [min([nopeaks;-eps(1)]),max([nopeaks;eps(1)])];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_tanacceleration.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_tanacceleration.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_tanacceleration.pdf'),'-fillpage')
  end
 end
 
 if ismember(5,plotfig)
  figure(5)
  hold on
  plot(lambda(2:end),An(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('normal acceleration');
  title('normal acceleration');
  grid on
  %nopeaks = (An(An<10*mean(An(:))));
  %bbb = gca();
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_noracceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_noracceleration.png'))
  end
 end
 
 if ismember(6,plotfig)
  figure(6)
  hold on
  plot(lambda(2:end),A0(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('total acceleration');
  title('total acceleration');
  grid on
  %nopeaks = (A0(A0<10*mean(A0(:))));
  %bbb = gca();
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_totacceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_totacceleration.png'))
  end
 end
 
 if ismember(7,plotfig)
  figure(7)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %plot(lambda(4:end),abs(TAU(4:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  plot(lambda(1:end),(res(k3).TAU(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o',2'Color',colo);
  xlabel('lambda');
  ylabel('second Frenet-curvature $\kappa_2$','Interpreter','latex');
  title(modelfilename,'Interpreter','none');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %nopeaks = (TAU(TAU<10*mean(TAU(:))));
  bbb = gca();
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  %line(xlim(), [1,1], 'LineWidth', eps(0), 'Color', 'k');
  xl = xlim; % Get current limits.
  xlim([0, xl(2)]); % Replace lower limit only with a y of 0.
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_torque.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_torque.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_torque.pdf'),'-fillpage')
  end
  disp(median(abs(res(k3).TAU(2:end)),'omitnan'));
 end
 
 if ismember(8,plotfig)
  figure(8)
  hold on
  plot(lambda(2:end),(sinpsi(2:end)),'LineStyle','-','Marker','x','LineWidth',1.5); %,'Color',col
  plot(lambda(2:end),(cosmu(2:end)),'LineStyle','-','Marker','o','LineWidth',1.5); %,'Color',col
  xlabel('lambda');
  ylabel('absolute values of: sin(psi) and cos(mu)');
  title('absolute values of: sin(psi) and cos(mu)');
  grid on
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_SinCos.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_SinCos.png'))
  print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_SinCos.pdf'),'-fillpage')
  end
 end
 
 if ismember(8.1,plotfig)
  figure(81)
  hold on
  plot(lambda(2:end),cosmu(2:end)./(sinpsi(2:end)),'LineStyle','-','Marker','x','LineWidth',1.5); %,'Color',col
  xlabel('lambda');
  ylabel('cosmu(2:end)./(sinpsi(2:end)');
  title('cosmu(2:end)./(sinpsi(2:end)');
  grid on
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_Sin-Cos.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_Sin-Cos.png'))
  print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_Sin-Cos.pdf'),'-fillpage')
  end
 end
 
 if ismember(11,plotfig) && numel(LAM)>0
  figure(11);
  hold on
  %        plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none');%,'Color',colo);
  %        plot(lambda(2:end), LAM(2:end),'LineStyle','--','Marker','none');%,'Color',colo);
  aaa = gca();
  %plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  plot(lambda(2:end), LAM(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  if size(aaa.Children,1)<2
   %plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
   %plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
  end
  xlabel('Lambda');
  %        ylabel('Lambda + chi | or | chi');
  ylabel('LAM=chi-lambda');
  grid on
  %        legend('chi','eigenvalue','lambda','zero')
  title(modelfilename)
  ylim([min(0,real(LAM(end))),max(min(real(LAM(end)),2),1)])
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_LAM.pdf'),'-dpdf')
  end
  if ismember(11.1,plotfig)
   ylim([-.1,.1])
   if model.savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAMX.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAMX.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAMX.pdf'),'-fillpage')
   end
  end
 end
 if ismember(-11,plotfig)
  if max(LAM(2:end))<0.01
   disp(LAM(2:end))
  else
   disp([lambda(2:end), LAM(2:end)])
  end
 end
 
 
 
 if ismember(12,plotfig)
  figure(12);
  hold on
  aaa = gca();
  plot(lambda(2:end),real(LAM(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('real(EW)=real(chi-lambda)');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  aaa.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM12.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM12.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM12.pdf'),'-fillpage')
  end
 end
 
 if ismember(13,plotfig)
  figure(13);
  hold on
  %        plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none');%,'Color',colo);
  %        plot(lambda(2:end), LAM(2:end),'LineStyle','--','Marker','none');%,'Color',colo);
  aaa = gca();
  if size(aaa.Children,1)<2
   %plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
   %plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
  end
  %plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  plot(lambda(2:end),abs(LAM(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  %        ylabel('Lambda + chi | or | chi');
  ylabel('abs(EW)=abs(chi-lambda)');
  title(modelfilename)
  grid on
  %        legend('chi','eigenvalue','lambda','zero')
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM13.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM13.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM13.pdf'),'-fillpage')
  end
 end
 
 
 
 
 if ismember(14,plotfig)
  figure(14);
  %set(gca, 'YScale', 'log')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   %grid minor
  end
%   xlabelJK='$\lambda$';
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='radius of the first Frenet-curvature $\rho$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
%   plot(lambda(2:end),(res(k3).RHO2(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(xValload,(res(k3).RHO2(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho14.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho14.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho14.pdf'),'-dpdf')
   %plotitJK(lambda(2:end),res(k3).RHO2(2:end),'Output/Figures/',xlabelJK,ylabelJK,strcat(modelfilename,'_rho14'),NaN,NaN,cfig) 
  end
  disp(median(res(k3).RHO2(2:lamlast),'omitnan'));
 end
 if ismember(-14,plotfig)
  disp([lambda(2:end),res(k3).RHO2(2:end)])
 end
 

 if ismember(15,plotfig)
  figure(15);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  x=xValfullload0; %model.fulllambda(1:end)
  y4=real(model.fullEV(k3,1:numel(x)));
%   plot(model.fulllambda(1:end),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  plot(x,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if k3==resEWs(1) && main.colorshift==0
   plot(x,zeros(size(x)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
%   xlabel('Lambda');
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='eigenvalue';
  ylabel(yLabel);
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if strcmp(main.typeofanalysis,'KNL2')
   maxy4=max(y4);
    bbb = gca();
    obergrenze=min(2,max(round(abs(maxy4),2,'significant'),maxy4));
    bbb.YLim = [-.2,obergrenze];
    cfig('yMax')=obergrenze;
  end
  if main.savefigures>=true
   figname='_LAM15';
   dianame=strcat(modelfilename,figname);
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   %    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM15.svg'))
   %    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM15.png'))
   if main.savefigures==2
    remove(cfig,{'xMin','Cxticks'});
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('yMin')=0;
    cfig('LineStyle')=MyLines{main.colorshift+1};
    plotitJK(x,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8015);
   end
  end
 end
 
 
 if ismember(16,plotfig)
  figure(16);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  
  plot(lambda(3:end),(res(k3).EWd2l(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\frac{(EW(x-\Delta)-2EW(x)+EW(x+\Delta)}{(\Delta^2)}$','Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
%   aaa = gca(); 
  %aaa.XLim = [0 inf];
  %set(aaa, 'YScale', 'log');aaa.YLim = [1e-4 1e1];
  %aaa.YLim = [-inf 1];
  %
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_EWd2l.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_EWd2l.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_EWd2l.pdf'),'-fillpage')
  end
%   if ismember(16.1,plotfig)
%    ylim([-inf 0])
%    if model.savefigures==true
%     print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_EWd2lneg.svg'))
%     print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_EWd2lneg.png'))
%     print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_EWd2lneg.pdf'),'-fillpage')
%    end
%   end
 end
 
 if ismember(17,plotfig)
  figure(17);
  hold on
  
  plot(lambda(3:end),(res(k3).OC6(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('OC6','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_OC6.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_OC6.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_OC6.pdf'),'-fillpage')
  end
 end
 
 if ismember(18,plotfig)
  figure(18);
  hold on
  
  plot(lambda(2:end),(res(k3).OC7(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('2*r(x)*[r(x+$\Delta$)+r(x-$\Delta$)]-r(x+$\Delta$)*r(x-$\Delta$)-3','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_OC7.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_OC7.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_OC7.pdf'),'-fillpage')
  end
 end
 
 if ismember(19,plotfig)
  figure(19);
  hold on
  aaa = gca();
  plot(lambda(2:end),imag(LAM(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('imag(EW)');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  aaa.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM12.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM12.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM12.pdf'),'-fillpage')
  end
 end
 
 if ismember(20,plotfig)
  figure(20);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(3:end-1),res(k3).rhoDach(3:end-1),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda $\lambda$');
  ylabel('$\hat{\rho}$ [-]','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rhorho3D.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rhorho3D.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_rhorho3D.pdf'),'-fillpage')
  end
 end
 
 if ismember(21,plotfig)
  figure(21);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  plot(lambda(2:end),res(k3).SinGamma(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\sin(\gamma)$','Interpreter','latex');
  %bbb = gca();
  %bbb.YLim = [0.0,1];
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_CosPhi.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_CosPhi.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_CosPhi.pdf'),'-fillpage')
  end
 
 end
  if ismember(211,plotfig)
   figure(211);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[0   0   XBreite   YHohe]); 
   hold on
   plot(lambda(2:end),res(k3).X1(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x1$','Interpreter','latex');
   ylabel('$\mathbf{r}_1\cdot\mathbf{r}_1"+1$','Interpreter','latex');
   %title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if model.savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X1.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X1.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X1.pdf'),'-fillpage')
   end
   figure(212);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[XBreite   0   XBreite   YHohe]); 
   hold on
   plot(lambda(2:end),res(k3).X2(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x2$','Interpreter','latex');
   title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if model.savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X2.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X2.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X2.pdf'),'-fillpage')
   end
   figure(213);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[2*XBreite   0   XBreite   YHohe]); 
   hold on
   plot(lambda(2:end),res(k3).X3(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('x3','Interpreter','latex');
   title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if model.savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X3.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X3.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X3.pdf'),'-fillpage')
   end
   figure(214);
   set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'XAxisLocation','origin')
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[3*XBreite   0   XBreite   YHohe]); 
   hold on
   set(gca, 'YScale', 'log')
   plot(lambda(2:end),res(k3).X4(2:end),'LineStyle','-','Marker','o','LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x_4$','Interpreter','latex');
   ylabel('$\mathbf{r}_1\cdot\mathbf{r}_1^{,,}+1$','Interpreter','latex');
   %title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if model.savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X4.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X4.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X4.pdf'),'-fillpage')
   end
  end
 
 if ismember(22,plotfig)
  figure(22);
  hold on
  
  plot(lambda(3:end),(res(k3).OC5(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('OC5','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_OC5.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_OC5.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_OC5.pdf'),'-fillpage')
  end
 end
 
 if ismember(23,plotfig)
  figure(23);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  xlabel('lambda');
  ylabel('$Hypothese=rho^2 sqrt(1+tau^2)$','Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0.0,1];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(lambda(1:end),(res(k3).HYPO(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_HYPO.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_HYPO.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_HYPO.pdf'),'-dpdf')
  end
 end
  
 
 
 if ismember(24,plotfig)
  figure(24);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  xlabel('lambda');
  ylabel('$\tau$','Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0.0,10];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(lambda(1:end),(res(k3).TAU(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_tau.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_tau.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_tau.pdf'),'-dpdf')
  end
 end
 if ismember(24.3,plotfig)
  figure(243);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==24.3)   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  xlabel('lambda');
  ylabel('$\tau$','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [-10,10];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(lambda(1:end),(res(k3).TauDach(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_tau_3D.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_tau_3D.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_tau_3D.pdf'),'-dpdf')
  end
 end
 
 if ismember(25,plotfig)
  figure(25);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  xlabel('lambda');
  ylabel('$\cos(\gamma)$','Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0.0,10];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(1:end),(res(k3).TAU(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color','red');
  %plot(lambda(1:end),(res(k3).ABSrddd(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(lambda(1:end),(res(k3).CosGamma(1:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  
  %yticks(0:.05:1)
  bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_ABSrddd.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_ABSrddd.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_ABSrddd.pdf'),'-dpdf')
  end
 end
 
 if ismember(26,plotfig)
  figure(26);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(2:end),res(k3).RXB(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('Skalarpodukt von $\mathbf{r}$ mit $\mathbf{b}$: ($\mathbf{r}\cdot\mathbf{b})$','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_RxB_26.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_RxB_26.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_RxB_26.pdf'),'-fillpage')
  end
  disp(median(res(k3).RXB(2:end),'omitnan'));
 end
 
 if ismember(27,plotfig)
  figure(27);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(2:end),res(k3).RXBDach(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$dot(\mathbf{\hat{r}},\mathbf{\hat{b}})$','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_RxB3D_27.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_RxB3D_27.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_RxB3D_27.pdf'),'-fillpage')
  end
 end
 
 if ismember(28,plotfig)
  figure(28);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(2:end),abs(res(k3).Rconst(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$|\mathbf{r}(0)\cdot\mathbf{r}(\lambda)|$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'Rconst.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'Rconst.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'Rconst.pdf'),'-fillpage')
  end
 end
 
 if ismember(29,plotfig)
  figure(29);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(2:end),res(k3).EBENE(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$EBENE$','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'EBENE.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'EBENE.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'EBENE.pdf'),'-fillpage')
  end
 end
 
 if ismember(30,plotfig)
  figure(30);
  %set(gca, 'YScale', 'log')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   %grid on
   %grid minor
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='radius of curvature $\rho_1$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0.001,1];
  %bbb.YLim = [0,1];
  title(modelfilename,'Interpreter','none')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(xValload,(res(k3).RHO2(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 60];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho30.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho30.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho30.pdf'),'-dpdf')
   if export
    %cfig('ylog')=true;
    %cfig('yMin')=0.0001;
    %cfig('yMax')=1;
    %close all
    plotitJK(xValload,res(k3).RHO2(2:end),'Output/',xlabelload,ylabelJK,strcat(modelfilename,'_rho30'),cfig)% %#ok<UNRCH>
   end
  end
  %disp(median(res(k3).RHO2(2:end),'omitnan'));
 end
 
 if ismember(31,plotfig)
  figure(31);
  %set(gca, 'YScale', 'log')
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   %grid minor
  end
  %xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='1-$\rho$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0.001,1];
  title(modelfilename,'Interpreter','none')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %plot(xValload,1-(res(k3).RHO2(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho31.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho31.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho31.pdf'),'-dpdf')
   if export
    plotitJK(xValload,1-res(k3).RHO2(2:end),'Output/Figures/',xlabelload,ylabelJK,strcat(modelfilename,'_rho31'),NaN,NaN,cfig)% %#ok<UNRCH>
   end
  end
  %disp(median(res(k3).RHO2(2:end),'omitnan'));
 end
 
 if ismember(32,plotfig)
  figure(32);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(xValload,res(k3).RXB(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  %xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$(\mathbf{r}\cdot\mathbf{e}_3)$';
  ylabel(ylabelJK,'Interpreter','latex');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_RxB_32.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_RxB_32.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_RxB_32.pdf'),'-fillpage')
   if export
    cfig('yMin')=-1;% %#ok<UNRCH>
    %plotitJK(xValload,res(k3).RXB(2:end),'Output/Figures/',xlabelload,ylabelJK,strcat(modelfilename,'_RxB_32'),cfig)
   end
  end
 end
 
 if ismember(33,plotfig)
  figure(33);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  phimin=min(res(k3).PHIR(2:end));
  phimax=max(res(k3).PHIR(2:end));
  if phimin==phimax
   bbb = gca();
   bbb.YLim = [phimin-eps(1),phimax+eps(1)];
  end
  plot(lambda(2:end),res(k3).PHIR(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$arccos(|\mathbf{r}(0)\cdot\mathbf{r}(\lambda)|)$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'PHIR.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'PHIR.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'PHIR.pdf'),'-fillpage')
  end
 end
 
 if ismember(34,plotfig)
  figure(34);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  plot(lambda(2:end),res(k3).DrhopDs(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\frac{\mathrm d\dot{\rho}}{\mathrm d\dot{s}}=-\kappa_2\,(\mathbf{r}\cdot\mathbf{e}_3)$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   %print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'PHIR.svg'))
   %print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'PHIR.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'DrhopDs34.pdf'),'-fillpage')
  end
  disp(median(-res(k3).DrhopDs(2:end),'omitnan'));
  disp(median(abs(res(k3).DrhopDs(2:end)),'omitnan'));
 end
 
 if ismember(35,plotfig)
  figure(35);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   plot(xValload,zeros(size(xValload)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  xValload0=[0;xValload];
  y4=real(model.fullEV(k3,1:end-3));
  plot(xValload0,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  xlabel('Lambda');
  ylabel('Eigenwert');
  title(modelfilename,'Interpreter','none')
  if strcmp(main.typeofanalysis,'KNL2')
   maxy4=max(y4);
    bbb = gca();
    obergrenze=min(2,round(abs(maxy4),2,'significant'));
    bbb.YLim = [-.2,obergrenze];
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM15.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM15.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM15.pdf'),'-fillpage')
  end
 end

 if ismember(36,plotfig)
  lastValue=min(size(model.fullEV,2),numel(xValfullload0));
  y4=real(model.fullEV(k3,1:lastValue));
  yLabel='eigenvalue';
  figname='_LAM36';
  dianame=strcat(modelfilename,figname);
  set(0, 'DefaultFigureWindowState', 'minimized');
  figure(36);
  xPlot=xValfullload0; %maybe if chrased?
  [idxs,obergrenze]=markMins(xPlot,y4);
  hold on
  if main.savefigures~=2
   %curfig=gca();
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   if k3==resEWs(1) && main.colorshift==0
    plot(model.fulllambda(1:end),zeros(size(model.fulllambda(1:end))),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
    grid on
    grid minor
   end
   %xPlot=xValfullload0; %model.fulllambda(1:end)
   
   plot(xPlot,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
   xlabel(xlabelload,'interpreter','latex');
   ylabel(yLabel);
   title(modelfilename,'Interpreter','none')
   if strcmp(main.typeofanalysis,'KNL2')
    maxy4=max(y4);
    bbb = gca();
    obergrenze=min(obergrenze,min(2,round(abs(maxy4),2,'significant')));
    if min(y4)<0
     bbb.YLim = [max(-.2,-.1*abs(obergrenze)),obergrenze];
    else
     bbb.YLim = [0,obergrenze];
    end
   end
   %markMins(xValfullload0,y4,gca());
   if model.savefigures==true
    %print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
    %print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   end
  end
  if main.savefigures==2
   %remove(cfig,{'xMin','Cxticks','yMax'});
   %    cfig('NaNs')='keep';
   %    cfig('closing')=NaN;
   %    cfig('lineColor')=colJK;
   %    cfig('yMin')=0;
   %    cfig('LineStyle')=MyLines{main.colorshift+1};
   %figure(8015);
   %obergr=curfig.CurrentAxes.YLim(2);
   %xlabel(xlabelload);
   cfig('wofuer')='TeX';
   cfig('yMax')=2.5e-4;%obergrenze;
   plotitJK(xPlot,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8036);
   %obergrenze=min(max(curfig.CurrentAxes.YLim(2),max(y4)),round(10*abs(min([min(y4) maxy4])),2,'significant'))
   %gca().YLim = [0,obergrenze];
   markMins(xPlot,y4);
   %xlabel(xlabelload);
  end
 end

 if ismember(37,plotfig)
  figure(37);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   %plot(model.fulllambda(1:end),zeros(size(model.fulllambda(1:end))),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  ls=res(1).stability_limit(end);
  x=xValload0;
  xlam=model.fulllambda(1:NrLoad0);
  y4=transpose(real(model.fullEV(k3,1:NrLoad0)))-(ls-xlam(:))/ls;
%   plot(x,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  plot(x,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
%   xlabel('Lambda');
  xlabel(xlabelload,'Interpreter','latex');
  ylabel('EW-$\frac{\lambda_S-\lambda}{\lambda_S}$','Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM37.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM37.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM37.pdf'),'-fillpage')
  end
 end

 if ismember(38,plotfig)
  figure(38);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$cos(\varphi)$ [-]';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xValload,(res(k3).cosPhiMangA(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  plot(xValload,(res(k3).cosPhiMangB(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'_cosphi38');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 if ismember(39,plotfig)
  figure(39);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$(d2rds2*Kt0_0*rm)$ [-]';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xValload,(res(k3).ZaelerB(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'ZaelerB39');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 if ismember(40,plotfig)
  figure(40);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$norm(d2rds2)$ [-]';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xValload,(res(k3).NormB1(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'NormB1_40');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 if ismember(41,plotfig)
  figure(41);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$norm(Kt0_0*rm)$ [-]';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xValload,(res(k3).NormB2(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'NormB2_41');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 if ismember(42,plotfig)
  figure(42);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$norm(rm)$ [-]';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(xValload,(res(k3).NormR(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'NormR_42');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 end %for k3
 
  if all(isnan(lambda))
   lambda=model.lambda0;
  end
 
 %model.savefigures=main.savefigures;
 colJK=MyColours{mod(main.colorshift,19)+1};
 cfig('lineColor')=colJK;
 %lambda1 = res.lambdaorg;
 %lambda = res(k3).lambda;
 if numel(lambda)<=41
  markJK=MyMarker(mod(main.colorshift,2)+1);
 else
  markJK='none';
 end
 Nr=main.colorshift;
 lineStyleJK=MyLines{mod(Nr,4)+1};
 cfig('LineStyle')=lineStyleJK;

 
 if ismember(900,plotfig)
  figure(900);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  lastvalue=min(numel(model.DetKtx),numel(lambda));
  relDetKtx=model.DetKtx(1:lastvalue)./1;%model.DetKtx(1);
  %model.lambdainput
%   minDetKtx=min(relDetKtx);
%   idxs=islocalmin(relDetKtx,'FlatSelection', 'all');
%   maxDetKtx=max(relDetKtx);
  plot(lambda,zeros(size(lambda)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambda(1:lastvalue),relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  markMins(lambda,relDetKtx);
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\det(K_T)$','Interpreter','latex');%/\det(K_T)_0
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'DetKt900.pdf'),'-fillpage')
  end
 end
 
 if ismember(901,plotfig)
  figure(901);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  DetKtx=model.DetKtx;
  lambdaplot=model.lambdainput;
  lambdaplot(DetKtx<0)=[];
  DetKtx(DetKtx<0)=[];
  relDetKtx=DetKtx.^(1/model.N0);
  %relDetKtx=relDetKtx./relDetKtx(1);
  plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambdaplot,relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\sqrt[N]{\det(K_T)}$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'sqrtDetKt901.pdf'),'-fillpage')
  end
 end
 
 if ismember(906,plotfig) || ismember(902,plotfig)
  h906=figure(906);
  if ismember(906,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h906, 'Visible', 'off');
  end
  hold on
  %y=cell2mat(model.arclengths);
  y4=1./cell2mat(model.dxidl);
  lambdaplot=model.load0;%lambda; lambdaplot=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\lambda/d\xi(JK) [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
%   [miny4,idx]=min(y4);
  idxs=islocalmin(y4,'FlatSelection', 'all');
  plot([0 0],[0 0],'LineWidth',eps(0));
  plot(lambdaplot(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  lambdaplot(idxs)
  %if k3==resEWs(1) && main.colorshift==0
  grid on
  grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldi906.pdf'),'-fillpage')
  end
 end

 if ismember(902,plotfig)
  figure(902);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthuJK);
  lambdaplot=model.load0;%lambda;
  NrPlot=min(numel(y),find(~isnan(lambdaplot),1,'last'));
  plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  xPlot=lambdaplot(1:NrPlot);
  yPlot=y(1:NrPlot);
  plot(xPlot,yPlot,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
%   xLabel=xlabelload;%'$\lambda$';
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='arc-length $\xi [m]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures>=true
   dianame=strcat(modelfilename,'xi902');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    cfig('yMax')=3;
    plotitJK(xPlot,yPlot,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8902);
    plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   end
  end
 end
 
 if ismember(903,plotfig)
  figure(903);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengths);
  %y4=cell2mat(model.arclengthsJK);
  lambdaplot=model.lambdainput;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  %plot(lambdaplot,y4(:,5),'LineStyle',':','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\xi [-]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
  grid on
  grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dxi903.pdf'),'-fillpage')
  end
 end
 
 if ismember(904,plotfig)
  figure(904);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %y=cell2mat(model.arclengths);
  y4=cell2mat(model.darclengthsJK);
  lambdaplot=model.lambdainput;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot,y4(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\xi(JK) [\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
  grid on
  grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dxi904.pdf'),'-fillpage')
  end
 end

 if ismember(905,plotfig)
  figure(905);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %y=cell2mat(model.arclengths);
  y4=cell2mat(model.dxidl);
  lambdaplot=model.lambdainput;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\xi(JK)/d\lambda [\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
  grid on
  grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dxidl905.pdf'),'-fillpage')
  end
 end
 
 if ismember(907,plotfig)
  figure(907);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthuJK);
  lambdaplot=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
%   plot(lambdaplot,cell2mat(model.arclengthurHM),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
 % plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
%   plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(y,lambdaplot(1:numel(y)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  ylabel('$\lambda$','Interpreter','latex');
  xlabel('arc-length $\xi [m]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'xi907.pdf'),'-fillpage')
  end
 end

 if ismember(909,plotfig) || ismember(908,plotfig)
  h909=figure(909);
  if ismember(909,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h909, 'Visible', 'off');
  end
  hold on
  %y=cell2mat(model.arclengths);
  y4=1./cell2mat(model.dvdl);
  lambdaplot=model.load0;%=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\lambda/dv [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  idxs=markMins(lambdaplot,y4);
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw909.pdf'),'-fillpage')
  end
 end
 
 if ismember(908,plotfig)
  figure(908);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=cell2mat(model.vMaxJK);
  lambdaplot=model.load0;%=lambda;
  x=lambdaplot(1:numel(y4));
  plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  plot(x,y4,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel(xlabelload,'Interpreter','latex');%'$\lambda$'
  yLabel='$v [m]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures>=true
   dianame=strcat(modelfilename,'w908');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
%     remove(cfig,{'yMax'});
    plotitJK(x,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8908);
    plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   end
  end
 end

 if ismember(910,plotfig)
  figure(910);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  lastvalue=min(numel(model.DetKtx),numel(lambda));
  relDetKtx=model.DetKtx(1:lastvalue)./1;%model.DetKtx(1);
  %model.lambdainput
%   minDetKtx=min(relDetKtx);
  idxs=islocalmin(relDetKtx,'FlatSelection', 'all');
%   maxDetKtx=max(relDetKtx);
  plot(lambda,zeros(size(lambda)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambda(1:lastvalue),relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  plot(lambda(idxs),relDetKtx(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
%   lambda(idxs)
%   minDetKtxs=relDetKtx(idxs);
%   minDetKtx=min([minDetKtxs(minDetKtxs>0);maxDetKtx]);
%   if minDetKtx<.01*maxDetKtx || any(minDetKtxs<0)
%    untergrenze=round((minDetKtx),2,'significant');
%    bbb = gca();
%    obergrenze=min(bbb.YLim(2),100*abs(untergrenze));
%    bbb.YLim = [0,obergrenze];
%   end
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\det(K_T)$','Interpreter','latex');%/\det(K_T)_0
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'DetKt900.pdf'),'-fillpage')
  end
 end

 if ismember(911,plotfig)
  figure(911);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.wMaxJK);
  lambdaplot=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
%   plot(lambdaplot,cell2mat(model.arclengthurHM),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
 % plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
%   plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(y)),y,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$w [m]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'w911.pdf'),'-fillpage')
  end
 end

 if ismember(912,plotfig)
  figure(912);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %y=cell2mat(model.arclengths);
  y4=1./cell2mat(model.dwdl);
  lambdaplot=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\lambda/dw [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
%   [miny4,idx]=min(y4);
  idxs=islocalmin(y4,'FlatSelection', 'all');
  plot([0 0],[0 0],'LineWidth',eps(0));
  plot(lambdaplot(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  lambdaplot(idxs)
  %if k3==resEWs(1) && main.colorshift==0
  grid on
  grid minor
  %end
  y4s=y4(idxs);
  maxy4=max(y4);
  miny4=min([min(y4s(y4s>0)) maxy4]);
  if miny4<.01*maxy4 || any(miny4<0)
   untergrenze=(miny4);
   bbb = gca();
   obergrenze=min(bbb.YLim(2),round(100*abs(untergrenze),2,'significant'));
   bbb.YLim = [0,obergrenze];
  end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw909.pdf'),'-fillpage')
  end
 end

 if ismember(913,plotfig)
  lastvalue=min(numel(model.DetKtx),numel(xValload0));
  relDetKtx=model.DetKtx(1:lastvalue)./model.DetKtx(1);
  xPlot=xValload0;%lambda
  yPlot=relDetKtx;
  yLabel='$\det(K_T)/\det(K_T)_0$';
  if main.savefigures==2
  else
   figure(913);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   hold on
   plot(xPlot,zeros(size(xPlot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   plot(xPlot(1:lastvalue),yPlot,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   markMins(xPlot,yPlot);
   xlabel(xlabelload,'Interpreter','latex');
   ylabel(yLabel,'Interpreter','latex');%
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   bbb = gca();
   if any(min(relDetKtx)<-1)
    bbb.YLim = [-1,bbb.YLim(2)];
   end
   if any(max(relDetKtx)>2)
    bbb.YLim = [bbb.YLim(1),2];
   end
  end
  dianame=strcat(modelfilename,'DetKt913');
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  elseif main.savefigures==2
   cfig('yMax')=2;
   figure(8913)
   hold on 
   markMins(xPlot,yPlot);
   plotitJK(xPlot,yPlot,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8913);
   remove(cfig,{'yMax'});
   %markMins(xPlot,yPlot,gca());
  end
 end
 
 if ismember(914,plotfig)
  figure(914);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  lastvalue=min(numel(model.DetKtx),numel(xValload0));
  maxDet=max(model.DetKtx(1:lastvalue));
  relDetKtx=model.DetKtx(1:lastvalue)./maxDet;
  %xplot=lambda(1:lastvalue);
%   xplot=xValload0;
  y4=relDetKtx;
  plot(lambda,zeros(size(lambda)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(xValload0,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
%   xLabel='$\lambda$';
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='$\frac{\det(K_T)}{\max(\det(K_T))}$';
  ylabel(yLabel,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if main.colorshift==0
   grid on
   grid minor
  end
  bbb = gca();
  if any(min(relDetKtx)<-.8)
   bbb.YLim = [-.8,bbb.YLim(2)];
  end
  if any(max(relDetKtx)>2)
   bbb.YLim = [bbb.YLim(1),2];
  end
  if main.savefigures==true
   figname='DetKt914';
   dianame=strcat(modelfilename,figname);
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    remove(cfig,{'xMin'});
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('Faktor')=1;
    %    cfig('yMin')=-.2;
    plotitJK(xValload0,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8914);
   end
  end
 end
 
 if ismember(915,plotfig)
  h915=figure(915);
  if ismember(915,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h915, 'Visible', 'off');
  end
  hold on
  %y=cell2mat(model.arclengths);
  y4=1./cell2mat(model.dvdl);
  lambdaplot=model.load0;%=lambda;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthsurHM','arclengthsuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\lambda/dv [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  markMins(lambdaplot,y4,gca());
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw915.pdf'),'-fillpage')
  end
 end
 
 if ismember(916,plotfig)
  figure(916);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=cell2mat(model.uMaxJK);
  lambdaplot=model.load0;%=lambda
  x=lambdaplot(1:numel(y4));
%   plot(lambdaplot(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  plot(x,y4,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='$u [m]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
   grid on
   grid minor
  if main.savefigures>=true
   dianame=strcat(modelfilename,'w916');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    plotitJK(x,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,8916);
%     plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   end
  end
 end


end %function
