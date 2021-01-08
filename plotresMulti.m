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
 % 24 tau
 % 25 CosGamma
 % 26 r*b
 % 27 r*b (3D)
 
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
  modelfilename=model.filename;
 else
  modelfilename='unknown';
 end
 if ~exist('plotfig','var')
  plotfig=[1,2];
 end
 %MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};

 %maxlam = max(lambda);
 %maxlam = 1.0;
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
 model.savefigures=false;
 for k3=resEWs
  %k3 %#ok<NOPRT>
  if k3==resEWs(end)
   model.savefigures=main.savefigures;
  end
  colJK=MyColours{mod(k3+main.colorshift-1,19)+1};
 %lambda1 = res.lambdaorg;
 lambda = res(k3).lambda;
 if numel(lambda)<=41
  markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
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
 TAU = res(k3).TAU;
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
 

 

 
 lambda = reshape(lambda,length(lambda),1);
 LAM = reshape(LAM,length(LAM),1);
 %LAM = LAM/maxlam;
 if ismember(1,plotfig)
  figure(1);
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
  ylabel('chi=real(EW)+lambda');
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
  bbb.YLim = [0.0,1.1];
  bbb.XLim = [0 inf];
  %col = bbb.Children(1).Color;
  
  %        plot(lambda(kl),RHO(kl),'LineStyle','none','Marker','x','LineWidth',1.5,'Color',col);
  grid on
  title(modelfilename)
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
  plot(lambda(4:end),abs(TAU(4:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  %       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o',2'Color',colo);
  xlabel('lambda');
  ylabel('torque');
  title(modelfilename);
  grid on
  %nopeaks = (TAU(TAU<10*mean(TAU(:))));
  %bbb = gca();
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  %line(xlim(), [1,1], 'LineWidth', eps(0), 'Color', 'k');
  xl = xlim; % Get current limits.
  xlim([0, xl(2)]); % Replace lower limit only with a y of 0.
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_torque.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_torque.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_torque.pdf'),'-fillpage')
  end
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
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  xlabel('lambda');
  ylabel('$\rho$','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0.0,1];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(lambda(3:end),(res(k3).RHO2(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  bbb.XLim = [0 inf];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho14.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho14.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho14.pdf'),'-dpdf')
  end
 end
 if ismember(-14,plotfig)
  disp([lambda(2:end),res(k3).RHO2(2:end)])
 end
 
 
 if ismember(15,plotfig)
  figure(15);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %aaa = gca();
  plot(model.fulllambda(1:end),real(model.fullEV(k3,1:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  xlabel('Lambda');
  ylabel('Eigenwert');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %aaa.XLim = [0 inf];
  %aaa.YLim = [-1 1];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM15.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM15.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM15.pdf'),'-fillpage')
  end
 end
 
 
 if ismember(16,plotfig)
  figure(16);
  hold on
  
  plot(lambda(3:end),(res(k3).EWd2l(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\frac{(EW(x-\Delta)-2EW(x)+EW(x+\Delta)}{(\Delta^2)}$','Interpreter','latex');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %aaa = gca(); %#ok<NASGU>
  %aaa.XLim = [0 inf];
  %set(aaa, 'YScale', 'log');aaa.YLim = [1e-4 1e1];
  %aaa.YLim = [-0.04 0.04];
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
  title(modelfilename)
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
   title(modelfilename)
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
   ylabel('$dot(\mathbf{r},\mathbf{b})$','Interpreter','latex');
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
   plot(lambda(2:end),res(k3).X4(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$\frac{\partial \dot{\rho}}{\partial \dot{s}}=-\tau\,\left(\mathbf{r}\cdot\mathbf{b}\right)$','Interpreter','latex');
   title(modelfilename)
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
  bbb.YLim = [0.0,10];
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
  ylabel('$dot(\mathbf{r},\mathbf{b})$','Interpreter','latex');
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
 
 end %for k3
end %function
