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
<<<<<<< HEAD
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
 
 
 %900 detKt
 
 
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
   cfig = containers.Map;
   cfig('ylog')=false;
   cfig('yMin')=0;
   cfig('yMax')=1;
   cfig('xMin')=0;
   cfig('line')='-';
   cfig('wofuer')='TeX';
   cfig('lineColor')=[0, 0.4470, 0.7410];
   cfig('MyMarker')='none';
   cfig('LineWidth')=3;
   cfig('FontSize')=17.5;
   cfig('order')=NaN;
   cfig('closing')=NaN;
 if export && isfield(model,'load')
  lengthlam=min(lengthlam,numel(model.load)+1);
  if strcmp(model.filename(1),'p') %%#ok<UNRCH>
   xlabelload='bending moment $M$ [kN\,m]';
   xValload=model.load(1:lengthlam-1)/1000;
   cfig('xMax')=5000;
  elseif strcmp(model.filename(1),'T')
   xlabelload='line load $p$ [kN/cm]';
   xValload=model.load(1:lengthlam-1)/100000;
   cfig('xMax')=70;
   cfig('Cxticks')=0:10:70;
  elseif strcmp(model.filename(1),'e')
   xlabelload='normal force $N$ [kN]';
   xValload=model.load(1:lengthlam-1)/1000;
   cfig('xMax')=1600;
   cfig('Cxticks')=0:200:1600;
  else
   xlabelload='load';
   xValload=model.load(1:lengthlam-1);
  end
 end
 
=======
 % 24 tau
 % 25 CosGamma
 % 26 r*b
 % 27 r*b (3D)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
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
<<<<<<< HEAD
  if isfield(model,'filename')
   modelfilename=model.filename;
  else
   modelfilename=model.TestCase;
  end
=======
  modelfilename=model.filename;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 else
  modelfilename='unknown';
 end
 if ~exist('plotfig','var')
  plotfig=[1,2];
 end
 %MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};

 %maxlam = max(lambda);
 %maxlam = 1.0;
<<<<<<< HEAD
 if ~exist('main','var')
  main.savefigures=false;
  main.colorshift=0;
 end
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
 model.savefigures=false;%might get overwritten for k3==resEWs(end)
=======
 model.savefigures=false;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 for k3=resEWs
  %k3 %#ok<NOPRT>
  if k3==resEWs(end)
   model.savefigures=main.savefigures;
  end
  colJK=MyColours{mod(k3+main.colorshift-1,19)+1};
<<<<<<< HEAD
  %lambda1 = res.lambdaorg;
  lambda = res(k3).lambda;
  if numel(lambda)<=41
   markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
  else
   markJK='none';
  end
=======
 %lambda1 = res.lambdaorg;
 lambda = res(k3).lambda;
 if numel(lambda)<=41
  markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
 else
  markJK='none';
 end
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
% markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
 %V = res.V;
 %A = res.A;
 S = res(k3).S;
 A0 = res(k3).A0;%...total accerlation
 At = res(k3).At;
 An = res(k3).An;
 RHO = res(k3).RHO;
 %RHO2 = res.RHO2;
<<<<<<< HEAD
 %TAU = res(k3).TAU;
=======
 TAU = res(k3).TAU;
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
 

 

 
<<<<<<< HEAD
 lambda = reshape(lambda,numel(lambda),1);
=======
 lambda = reshape(lambda,length(lambda),1);
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 LAM = reshape(LAM,length(LAM),1);
 %LAM = LAM/maxlam;
 if ismember(1,plotfig)
  figure(1);
<<<<<<< HEAD
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  ylabel('EW');
=======
  ylabel('chi=real(EW)+lambda');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  %bbb.YLim = [0.0,1.1];
=======
  bbb.YLim = [0.0,1.1];
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  bbb.XLim = [0 inf];
  %col = bbb.Children(1).Color;
  
  %        plot(lambda(kl),RHO(kl),'LineStyle','none','Marker','x','LineWidth',1.5,'Color',col);
  grid on
<<<<<<< HEAD
  title(modelfilename,'Interpreter','none')
=======
  title(modelfilename)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
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
=======
  plot(lambda(4:end),abs(TAU(4:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  %       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o',2'Color',colo);
  xlabel('lambda');
  ylabel('torque');
  title(modelfilename);
  grid on
  %nopeaks = (TAU(TAU<10*mean(TAU(:))));
  %bbb = gca();
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  %bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  %line(xlim(), [1,1], 'LineWidth', eps(0), 'Color', 'k');
  xl = xlim; % Get current limits.
  xlim([0, xl(2)]); % Replace lower limit only with a y of 0.
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_torque.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_torque.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_torque.pdf'),'-fillpage')
  end
<<<<<<< HEAD
  disp(median(abs(res(k3).TAU(2:end)),'omitnan'));
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  %set(gca, 'YScale', 'log')
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
<<<<<<< HEAD
   %grid minor
  end
  xlabelJK='$\lambda$';
  xlabel(xlabelJK,'Interpreter','latex');
  ylabelJK='radius of the first Frenet-curvature $\rho$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  title(modelfilename,'Interpreter','none')
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(lambda(2:end),(res(k3).RHO2(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 inf];
=======
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
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho14.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho14.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho14.pdf'),'-dpdf')
<<<<<<< HEAD
   %plotitJK(lambda(2:end),res(k3).RHO2(2:end),'Output/Figures/',xlabelJK,ylabelJK,strcat(modelfilename,'_rho14'),NaN,NaN,cfig) 
  end
  disp(median(res(k3).RHO2(2:lamlast),'omitnan'));
=======
  end
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 end
 if ismember(-14,plotfig)
  disp([lambda(2:end),res(k3).RHO2(2:end)])
 end
 
 
 if ismember(15,plotfig)
  figure(15);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
<<<<<<< HEAD
  plot(model.fulllambda(1:end),real(model.fullEV(k3,1:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  xlabel('Lambda');
  ylabel('Eigenwert');
  title(modelfilename,'Interpreter','none')
=======
  %aaa = gca();
  plot(model.fulllambda(1:end),real(model.fullEV(k3,1:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  xlabel('Lambda');
  ylabel('Eigenwert');
  title(modelfilename)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
<<<<<<< HEAD
  aaa = gca();
  %aaa.XLim = [0 inf];
  if strcmp(main.typeofanalysis,'KNL2')
   aaa.YLim = [-.2 1.2];
  end
=======
  %aaa.XLim = [0 inf];
  %aaa.YLim = [-1 1];
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM15.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM15.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM15.pdf'),'-fillpage')
  end
 end
 
 
 if ismember(16,plotfig)
  figure(16);
<<<<<<< HEAD
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
=======
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  hold on
  
  plot(lambda(3:end),(res(k3).EWd2l(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\frac{(EW(x-\Delta)-2EW(x)+EW(x+\Delta)}{(\Delta^2)}$','Interpreter','latex');
<<<<<<< HEAD
  title(modelfilename,'Interpreter','none')
=======
  title(modelfilename)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
<<<<<<< HEAD
  aaa = gca(); 
  %aaa.XLim = [0 inf];
  %set(aaa, 'YScale', 'log');aaa.YLim = [1e-4 1e1];
  aaa.YLim = [-inf 1];
=======
  %aaa = gca(); %#ok<NASGU>
  %aaa.XLim = [0 inf];
  %set(aaa, 'YScale', 'log');aaa.YLim = [1e-4 1e1];
  %aaa.YLim = [-0.04 0.04];
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  title(modelfilename,'Interpreter','none')
=======
  title(modelfilename)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
   ylabel('x3','Interpreter','latex');
=======
   ylabel('$dot(\mathbf{r},\mathbf{b})$','Interpreter','latex');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
   ylabel('$x_4$','Interpreter','latex');
=======
   ylabel('$\frac{\partial \dot{\rho}}{\partial \dot{s}}=-\tau\,\left(\mathbf{r}\cdot\mathbf{b}\right)$','Interpreter','latex');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  %bbb.YLim = [0.0,10];
=======
  bbb.YLim = [0.0,10];
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  ylabel('Skalarpodukt von $\mathbf{r}$ mit $\mathbf{b}$: ($\mathbf{r}\cdot\mathbf{b})$','Interpreter','latex');
=======
  ylabel('$dot(\mathbf{r},\mathbf{b})$','Interpreter','latex');
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
<<<<<<< HEAD
  disp(median(res(k3).RXB(2:end),'omitnan'));
 end
 
 if ismember(27,plotfig)
=======
 end
 
  if ismember(27,plotfig)
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
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
 
<<<<<<< HEAD
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
    plotitJK(xValload,res(k3).RHO2(2:end),'Output/Figures/',xlabelload,ylabelJK,strcat(modelfilename,'_rho30'),cfig)% %#ok<UNRCH>
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
 

 
 end %for k3
 

 %model.savefigures=main.savefigures;
 colJK=MyColours{mod(main.colorshift,19)+1};
 %lambda1 = res.lambdaorg;
 lambda = res(k3).lambda;
 if numel(lambda)<=41
  markJK=MyMarker(mod(main.colorshift,2)+1);
 else
  markJK='none';
 end

 
 if ismember(900,plotfig)
  figure(900);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  relDetKtx=model.DetKtx./1;%model.DetKtx(1);
  plot(model.lambdainput,relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\det(K_T)/\det(K_T)_0$','Interpreter','latex');
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
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'DetKt901.pdf'),'-fillpage')
  end
 end
 
=======
 end %for k3
>>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
end %function
