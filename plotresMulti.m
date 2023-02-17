function plotresMulti(res,model,plotfig,MyColours,MyMarker,resEWs,main)
%% Plot Graphs
%author:Johannes Kalliauer(©2020-2023) based on Michał Malendowski (©2019-2020)

%% Input
% res ... results from sortEigenValuesAndGetQuantities
% model ... reults from runEigenProblem
% plotfig ... which figures to plot
% MyColours ... which colours to use
% resEWs ... which EigenWalues
% main ... optional post-processing-parameters from the script


%% no Output (only matlab-figures, and png/pdf/svg/eps/wmf-files)

%% Structure
% * runEigenProblem ... run the Eingenvalue-Problem
%   * selectModel .. calls a function to create the input-file 
%   * AbaqusModelsGeneration.runAbaqus ... run the input-file in Abaqus
%   * AbaqusModelsGeneration.getStiffnessMatrices ... get the stiffness-matrix from Abaqus-results
%   * AbaqusModelsGeneration.getHistoryOutputFromDatFile ... get the nodal-results from Abaus
%   * runEigenProblemSub ... Run the core of the eigenvalue-Problem and saving it into model
%     * solveCLEforMinEigNew ... get one specific eigenvalue and the eigenvector
%   * runEigenProblemDispJK ... Posprocessing the displacements
% * sortEigenValuesAndGetQuantities ... does the calculation of \rho
% * plotresMulti ... Plots the requested graphs

%% recent-change
%2023-02-16 JK: if main.savefigures missing, setting default value to false

%% plotfig
 % 2...rho
 % 1...chi
 % 3..._velocity
 % 4_tanacceleration
 % 5_noracceleration
 % 6_totacceleration
 % 7_torque
 % 8_SinCos
 % 8.1 Sin/Cos
 % 11_ Eigenwert ylim([0,1]) alt LAM=chi-lambda
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
 %37 Eigenwert only nonlinearPart EW- (ls-l)/ls
 %38 cosPhiMang
 %39 ZaelerB
 %40 NennerB1
 %41 NennerB2
 %42 NormR
 %43 EVal LAM=chi-lambda
 %44 EWd1l
 %45 abs(model.fullEV)
 
 
 
 
 %900 detKt
 %901 sqrt[N]{detKt}
 %902 \xiJK
 %905 model.dxidl
 %906 dl/dxi
 %908 v [m]
 %909 d\lambda/dw
 %913 det/det0
 %914 det/max(det)
 %918 model.displacementsJK{step}(1,:)
 
 %943 model.rddotKtr
 %944 ZweirKtt
 %945 rdotKtr
 %947 res(k3).EWd2l(1:end)./model.rddotKtr(1:Xlength);
 %948 model.ZweirKtt(1:Xlength)./model.rddotKtr(1:Xlength)
 %949 model.Zwei_t_KB1_t(1:Xlength)./model.rddotKtr(1:Xlength);
 %951 rdotKtt
 %952 t_KB1_t
 %953 -2*model.t_KB1_t
 %958 model.rKt0r
 %972 model.rdotKtr / model.rKt0r
 
 %% Code
 
 xBezug=main.xBezug;
%n..normiert auf Reflast
% 1...ohne xValfulllambda0Mult;
% s...Stepnumber; 
% d..differenz: normiert, aber 0 bei Reflast; 
% alles andere... last
 
%  if numel(resEWs)>0 && ~strcmp(main.whichEV,'skip')
%   lambda = res(min(resEWs)).lambda;
%  else
  lambda=model.lambda;
%  end
 lengthlam=length(lambda);
 %lamlast=find(~isnan(lambda),1,'last');

 export=true;
 MyLines={'-','--','-.',':'};
   cfig = containers.Map;
   cfig('ylog')=false;
   cfig('yMin')=0;
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
   xValfulllambda0Mult=1;
   xlabelload='load';
   xValload=model.load(1:lengthlam-1);
   xValload0=model.load0;
   lamMin=NaN;
   %lamMax=NaN;
 if export && isfield(model,'load')
  lengthlam=min(lengthlam,numel(model.load)+1);
  if strcmp(model.filename(1),'p') %%#ok<UNRCH>%pureBending
   xlabelload='bending moment $M$'; % [kN\,m]
   xValload=model.load(1:lengthlam-1);%/1000;
   xValload0=model.load0;%/1000;
  elseif strcmp(model.filename(1),'T')%Stuetlinienbogen
   xlabelload='line load $p$ [kN/cm]';
   xValload=model.load(1:lengthlam-1);%/100000;
   xValload0=model.load0;%/100000;
   %cfig('Cxticks')=0:10:70;
   xValfulllambda0Mult=1/.333302;
  elseif strcmp(model.filename(2),'ec')%exzentrischer Druck alte Reflast
   xlabelload='normal force $N$ [N]';
   xValload=model.load(1:lengthlam-1);%/1000;
   xValload0=model.load0;%/1000;
   %cfig('Cxticks')=0:200:1600;
   xValfulllambda0Mult=1/2.2;%(500e3*(5)^2)/(pi^2*2.1*10^11*1.3198*10^-5);%1/2.2;
   lamMin=0;
   %lamMax=2;
  elseif strcmp(model.filename(2),'ex')%exzentrischer Druck neue Reflast
   xlabelload='normal force $N$ [N]';
   xValload=model.load(1:lengthlam-1);%/1000;
   xValload0=model.load0;%/1000;
   %cfig('Cxticks')=0:200:1600;
   xValfulllambda0Mult=1;%(1100e3*(5)^2)/(pi^2*2.1*10^11*1.3198*10^-5);
   lamMin=0;
   %lamMax=2;
  elseif strcmp(model.filename(1),'d')%IPE 400 bar subjected to a centric axial force at its midpoint
   if strcmp(model.filename(1:6),'detKtN')
    xValfulllambda0Mult=1;%https://www.wolframalpha.com/input/?i=2*%28pi%5E2*2.1*10%5E11*1.3198*10%5E-5%29%2F%285*.699156%29%5E2
   else
    xValfulllambda0Mult=1/66.545;%https://www.wolframalpha.com/input/?i=2*%28pi%5E2*2.1*10%5E11*2.18765*10%5E-4%29%2F%285*.699156%29%5E2
   end
   xlabelload='Force';
   xValload=model.load(1:lengthlam-1);
   xValload0=model.load0;
   lamMin=-1.5;
   %lamMax=1.5;
  elseif strcmp(model.filename(1),'B')%BendC...BendingCantilever
   xValfulllambda0Mult=1/0.57280;%/0.57280;
  elseif strcmp(model.filename(1),'c')%canti
   xlabelload='point load [N]';
  elseif strcmp(model.filename(1),'K')%Kreis_arch3D
   xValfulllambda0Mult=1/0.15;% rought estimate where the first complex eigenvalue gets negative
  else %everything unknown
   xlabelload='load';
   xValload=model.load(1:lengthlam-1);
   xValload0=model.load0;
   lamMin=NaN;
  end
%   xValload0=[0;xValload(1:end)];
  NrLoad0=numel(xValload0);
  xValfullload0=model.fullload0;
 end
 if sum(strcmp(fieldnames(model), 'xlabelloadname')) ~= 0
   xlabelload=model.xlabelloadname;
 end
 if ~exist('xlabelload','var')
  xlabelload='xlabelload';
 end
 if ~exist('xValload','var') || numel(xValload)==0
  warning('MyPrgm:Input','xValload not defined')
  xValload=model.fulllambda(2:end);
 end
 if ~exist('xValload0','var') || numel(xValload0)==0
  warning('MyPrgm:Input','xValload0 not defined')
  xValload0=[0;xValload];
 end

 
 if isempty(model.fulllambda)
  model.fulllambda=model.lambda0;
 end
 %LamdaObj=model.fulllambda(1:end)*xValfulllambda0Mult;
 lambdaLabel='$\lambda$';%lambdaLabel='load parameter $\lambda$';
 
 
  %% Figuredef
  FontSize0=(get(0,'defaultaxesfontsize')+get(0,'defaulttextfontsize'))/2;
  tmp=get(0,'ScreenSize');
  GetSceenX=tmp(3);
  %GetSceenY=tmp(4);
  %screenX= 1920;
  screenX=max(min(GetSceenX,3840),1024);
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
if isempty(resEWs)
 plotfig=plotfig(~((plotfig<=34).*(plotfig>0)));
end
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

 

 if ~exist('main','var')
  main.savefigures=false;
  main.colorshift=0;
 end
 if ~isfield(main,'savefigures'); main.savefigures=false; warning('MyPrgm:undefined','assuming not saving figures'); end
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
 %model.savefigures=false;%might get overwritten for k3==resEWs(end)
 %set(0, 'DefaultFigureWindowState', 'minimized');
 cfig('Faktor')=1;
 if lambda(1) == 0
  lambda0=zeros(numel(lambda),1);
  lambda0(1:end)=lambda;
 else
  if isnan(lambda(1))
   lambda0=NaN(numel(lambda)+1,1);
  else
   lambda0=zeros(numel(lambda)+1,1);
  end
  lambda0(2:end)=lambda;
 end
 NrFulllambda0=numel(model.fulllambda);

 if strcmp(xBezug,'n')%lambda nachträglich normalisiert
  xPlot0=lambda0(1:end)*xValfulllambda0Mult; %xValfullload0
  xPlot=model.fulllambda(2:end-3)*xValfulllambda0Mult; %xValfullload0
  FullxPlot0=model.fulllambda(1:end-3)*xValfulllambda0Mult;%different from xPlot0 if lambda(1)>>epsil
  myxlabelload='$\lambda$ gem Referenzlast';
  %xlim([.94 1.03])
 elseif strcmp(xBezug,'1') %Abaqus-Lambda
  xPlot0=model.lambdainput(1:end); %xValfullload0
  xPlot=xPlot0(2:end);%model.fulllambda(1:end); %xValfullload0
  FullxPlot0=model.fulllambda(1:end-3);%different from xPlot0 if lambda(1)>>epsil
  myxlabelload='$\lambda gem Abaqusinput$';
 elseif strcmp(xBezug,'s') %Lastscrhittnumber
  xPlot=1:NrFulllambda0; 
  xPlot0=[0 xPlot];
  myxlabelload='$Stepnumber$';
 elseif strcmp(xBezug,'d')%lambda differenz
  xPlot0=(lambda0(1:end)*xValfulllambda0Mult)-1; %xValfullload0
  xPlot=(model.fulllambda(1:end)*xValfulllambda0Mult)-1; %xValfullload0
  myxlabelload='$\frac{P-\bar{P}}{\bar{P}}$';
  %xlimJK=[-.1 .1];
  %xlimJK=[-1 1];
 elseif strcmp(xBezug,'P')
  xPlot0=xValload0(1:end);%x=lambda(3:end)
  xPlot=xValload0(2:end);%x=lambda(3:end)
  myxlabelload=xlabelload;
  FullxPlot0=model.fulllambda(1:end-3);%different from xPlot0 if lambda(1)>>epsil
 elseif strcmp(xBezug,'i') %Abaqus-Lambda
  xPlot0=model.lambdainput(1:end)/main.xBezugNr; %xValfullload0
  xPlot=xPlot0(2:end);%model.fulllambda(1:end); %xValfullload0
  myxlabelload='$\lambda gem Abaqusinput$';
 else
  error('MyPrgm:Unknown','xBezug not defined')
 end
 steps0=numel(xPlot0);
 if steps0-1<numel(xPlot)
  warning('MyPrgm:MissingValues','Not all requested lambdas have values try reruning modelprops.forceAbaqus=1')
  xPlot=xPlot(1:steps0-1);
 end
 xPlotLength=numel(xPlot);
 xPlot0Length=numel(xPlot0);
 tmpLast=min(xPlotLength,xPlot0Length-1);
  xPlotHalf=(xPlot0(1:tmpLast)+xPlot)/2;


 stepNR=numel(xPlot0);
 lambda = reshape(lambda,numel(lambda),1);
 for k3=resEWs
  %k3 %#ok<NOPRT>
  if k3==resEWs(end)
   savefigures=main.savefigures;
  else
   savefigures=false;
  end
  Nr=k3+main.colorshift-1;
  colJK=MyColours{mod(Nr,19)+1};
  cfig('lineColor')=colJK;
%   lineStyleJK=MyLines{mod(Nr,4)+1};
  %lambda1 = res.lambdaorg;

  if numel(lambda)<=41
   markJK=MyMarker(mod(Nr,2)+1);
  else
   markJK='none';
  end
% markJK=MyMarker(mod(k3+main.colorshift-1,2)+1);
  if  ~strcmp(main.whichEV,'skip')
   %lambdaRes = res(k3).lambda;
   S = res(k3).S;
   A0 = res(k3).A0;%...total accerlation
   At = res(k3).At;
   An = res(k3).An;
   RHO = res(k3).RHO;
   %TAU = res(k3).TAU;
   cosmu = res(k3).cosmu;
   sinpsi = res(k3).sinpsi;
   LAM = res(k3).LAM;
   LAM = reshape(LAM,length(LAM),1);
  end

 if ismember(1,plotfig)
  figure(1);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  aaa = gca();
  if size(aaa.Children,1)<2
   plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
   plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
  end
  plot(lambda(2:end), lambda(2:end)+real(LAM(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);% ,'Color',col %,'Color',colo);
  xlabel('Lambda');
  ylabel('EW');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  title(modelfilename)
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_chi.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_chi.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_chi.pdf'),'-dpdf')
  end
 end

%fignr=2; %fignr=02;
 if ismember(2,plotfig)
  figure(2);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]); 
  hold on
  Values=min(numel(lambda)-1,numel(res(k3).RHO)-1);
  plot(lambda(2:Values),RHO(2:Values),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  
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
  if savefigures==true
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
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_velocity.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_velocity.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_velocity.pdf'),'-fillpage')
  end
 end
 
 if ismember(3.3,plotfig)
  figure(330)
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
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_totacceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_totacceleration.png'))
  end
 end
 
 %plotfig=007;
 if ismember(7,plotfig)
  figure(7)
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  yData=res(k3).TAU(1:end);
  yNumber=numel(yData);
  %plot(lambda(4:end),abs(TAU(4:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5);%,'Color',col,'Color',colo);
  plot(lambda(1:yNumber),yData,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o',2'Color',colo);
  xlabel('lambda');
  ylabel('second Frenet-curvature $\kappa_2$','Interpreter','latex');
  title(modelfilename,'Interpreter','none');
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  bbb = gca();
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %line(xlim(), [1,1], 'LineWidth', eps(0), 'Color', 'k');
  xl = xlim; % Get current limits.
  xlim([0, xl(2)]); % Replace lower limit only with a y of 0.
  yl= ylim;
  if yl(2)>10
   ylim([0,min(10,yl(2))]);
  end
  if savefigures==true
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
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);

  hold on
  plot(lambda(2:end), LAM(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('LAM=chi-lambda');
  grid on
  title(modelfilename)
  ylim([min(-0.1,max(real(LAM(end)),-0.2)),max(min(real(LAM(end)),2),1)])
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
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);

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
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM12.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM12.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM12.pdf'),'-fillpage')
  end
 end
 
 if ismember(13,plotfig)
  figure(13);
  hold on
  plot(lambda(2:end),abs(LAM(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('abs(EW)=abs(chi-lambda)');
  title(modelfilename)
  grid on
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM13.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM13.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM13.pdf'),'-fillpage')
  end
 end
 
 
 
 %fignr=14;
 if ismember(14,plotfig)
  figure(14);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end  
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='radius of the first Frenet-curvature $\rho$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0,1];
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  if all(isnan(res(k3).RHO2)) % from sortEigenValuesAndGetQuantities.m
   warning('MyPrgm:Plot:Rho2:NaN','RHO2 is NaN')
  else
   Values=min(numel(xPlot)-1,numel(res(k3).RHO2)-1);
   plot(xPlot(2:Values),res(k3).RHO2(2:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  end
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho14.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho14.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho14.pdf'),'-dpdf')
  else
   title(modelfilename)
  end
  %disp(median(res(k3).RHO2(2:lamlast),'omitnan'));
 end

 if ismember(-14,plotfig)
  disp([lambda(2:end),res(k3).RHO2(2:end)])
 end
 
%fignr=15
 if ismember(15,plotfig)
  figure(15);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=real(model.fullEV(k3,:));
  stepNRtmp=min(stepNR,numel(y4));
  plot(xPlot0(1:stepNRtmp),y4(1:stepNRtmp),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   %plot(xPlot0(1:stepNRtmp),zeros(stepNRtmp),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  yLabel='$Re(\chi_1)$';
  ylabel(yLabel,'Interpreter','latex');
  xlabel(myxlabelload)
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if strcmp(main.typeofanalysis,'KNL2')
   maxy4=max([y4,1.0]);
    bbb = gca();
    obergrenze=min(2.5,max(round(abs(maxy4),2,'significant'),maxy4));
    bbb.YLim = [-1.2,obergrenze];
    cfig('yMax')=obergrenze;
  end
   figname='_LAM15';
   dianame=strcat(modelfilename,figname);
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   %print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM15.png'))
  elseif savefigures==2
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('LineStyle')=MyLines{main.colorshift+1};
    cfig('XAxisLocation')='bottom';
    cfig('xMin')=lamMin;
    cfig('yMax')=1;
    cfig('yMin')=0;
    h8015=figure(8015);
    plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8015);
    if k3==resEWs(1) && main.colorshift==0
    elseif k3==resEWs(1) && main.colorshift==1
     legend('B32OS','B32OSH')
     dianame=strcat(modelfilename,figname);
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))%
    end
  else
   title(modelfilename)
  end
 end
 
 %fignr=016;
 if ismember(16,plotfig)
  figure(16);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  Xlength=min(numel(xValload0),numel(res(k3).EWd2l));
  firstVal=5;
  y4=res(k3).EWd2l(firstVal:end);%real(model.fullEV(k3,1:numel(FullxPlot0)));%y4=real(model.fullEV(k3,1:end-3));
  resLamimag=imag(res(k3).LAM(firstVal:Xlength));
  if main.allowComplex==2
   y4(resLamimag==0)=NaN;
  elseif main.allowComplex==0
   y4(resLamimag~=0)=NaN;
  end
  plot(xPlot0(firstVal:Xlength),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  %xlabel('Lambda');  
  ylabel('$\frac{(EW(x-\Delta)-2EW(x)+EW(x+\Delta)}{(\Delta^2)}$','Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %ylim([-30,30])
  if model.savefigures ~= 2
   title(modelfilename)
  end
  if savefigures==true
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
 
 %fignr=19
 if ismember(19,plotfig)
  figure(19);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);

  hold on
  %xPlot=lambda0(1:end)*xValfulllambda0Mult;  
  y19=model.fullEV(k3,1:numel(FullxPlot0));
  y19imag=imag(y19);%y4=imag(LAM(2:end));
  yreal=real(y19);
  y19imag(isnan(yreal))=NaN;
  y4=y19imag;
  if main.allowComplex==2
   y4(y19imag==0)=NaN;
  end
  plot(FullxPlot0,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);

  xlabel('Lambda');
  ylabel('imag(EW)');
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  %aaa = gca();
  %aaa.XLim = [0 inf];
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM19.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM19.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM19.pdf'),'-fillpage')
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
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_CosPhi.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_CosPhi.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_CosPhi.pdf'),'-fillpage')
  end
 
 end
  if ismember(211,plotfig)
   figure(211);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[0   0   XBreite   YHohe]); 
   hold on
   LastVal=numel(res(k3).X1);
   plot(lambda(2:LastVal),res(k3).X1(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x1$','Interpreter','latex');
   ylabel('$\mathbf{r}_1\cdot\mathbf{r}_1"+1$','Interpreter','latex');
   %title(modelfilename)
   set(gca, 'YScale', 'log')
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X1.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X1.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X1.pdf'),'-fillpage')
   end
   %fignr=212
   figure(212);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[XBreite   0   XBreite   YHohe]); 
   hold on
   LastVal=numel(res(k3).X2);
   plot(lambda(2:LastVal),res(k3).X2(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x2$','Interpreter','latex');
   title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   ylim([-1.2,-0.8])
   if savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X2.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X2.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X2.pdf'),'-fillpage')
   end
   figure(213);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[2*XBreite   0   XBreite   YHohe]); 
   hold on
   LastVal=numel(res(k3).X3);
   plot(lambda(2:LastVal),res(k3).X3(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('x3','Interpreter','latex');
   title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X3.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X3.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X3.pdf'),'-fillpage')
   end
   figure(214);
   set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'XAxisLocation','origin')
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[3*XBreite   0   XBreite   YHohe]); 
   hold on
   set(gca, 'YScale', 'log')
   LastVal=numel(res(k3).X4);
   plot(lambda(2:LastVal),res(k3).X4(2:end),'LineStyle','-','Marker','o','LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('$x_4$','Interpreter','latex');
   ylabel('$\mathbf{r}_1\cdot\mathbf{r}_1^{,,}+1$','Interpreter','latex');
   %title(modelfilename)
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if savefigures==true
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
  ylabel('$Hypothese=rho (1+tau)$','Interpreter','latex');
  bbb = gca();
  bbb.YLim = [0.0,10];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  yData=res(k3).HYPO(1:end);
  yNumel=numel(yData);
  plot(lambda(1:yNumel),yData,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 inf];
  if savefigures==true
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
 
 %plotfig=26;
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
  if savefigures==true
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
  if savefigures==true
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
  tmplast=numel(res(k3).RHO2(2:end));
  %plot(lambda(3:end),res(k3).RHO2(3:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  plot(xValload(1:tmplast),(res(k3).RHO2(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  %bbb.XLim = [0 60];
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho30.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho30.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho30.pdf'),'-dpdf')
   if export && savefigures>1 
    h8030=figure(8030);
    plotitJK(xValload,res(k3).RHO2(2:end),'Output/',xlabelload,ylabelJK,strcat(modelfilename,'_rho30'),cfig,h8030)% %#ok<UNRCH>
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
  if savefigures==true
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
  %   if strcmp(xBezug,'n')%lambda
  %    xPlot=lambda(2:end)*xValfulllambda0Mult; %xValfullload0
  %    xlabel('$\lambda$','Interpreter','latex');
  %   elseif strcmp(xBezug,'1')
  %    xPlot=lambda(2:end); %xValfullload0
  %    xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  %   else
  %    xPlot=xValload0(3:Xlength);%x=lambda(3:end)
  %    xlabel(xlabelload,'Interpreter','latex');
  %   end
  tmplast=numel(res(k3).PHIR(2:end));
  plot(xPlot(1:tmplast),res(k3).PHIR(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  ylabel('$arccos(|\mathbf{r}(0)\cdot\mathbf{r}(\lambda)|)$','Interpreter','latex');
  %title(modelfilename,'Interpreter','none')
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'PHIR33.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'PHIR33.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'PHIR33.pdf'),'-fillpage')
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
 
 %plotfig=35;
 if ismember(35,plotfig)
  figure(35);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   %plot(xPlot0,zeros(size(xPlot0)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  y4=real(model.fullEV(k3,1:numel(FullxPlot0)));%y4=real(model.fullEV(k3,1:end-3));
   y19imag=imag(model.fullEV(k3,1:numel(xPlot0)));
  if main.allowComplex==2
   y4(y19imag==0)=NaN;
  elseif main.allowComplex==0
   y4(y19imag~=0)=NaN;
  end
  plot(FullxPlot0,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  %ylim([-.4 1])
 % xlim([-.8 -.4]) 
 %set(gca, 'YScale', 'log')
  xlabel(myxlabelload,'Interpreter','latex');
  ylabel('Re($\chi$)','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  if strcmp(main.typeofanalysis,'KNL2')
%    maxy4=max(y4);
%     bbb = gca();
%     obergrenze=max(1,min(700,round(abs(maxy4),2,'significant')));
%     untergr   =min(max(min(y4),-1),-1);
    %bbb.YLim = [untergr,obergrenze];%[.99,1];
  end
  %ylim([.1,.7])
  if savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM35.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM35.png'))
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_LAM35.pdf'),'-fillpage')
  end
 end

 if ismember(36,plotfig)
  lastValue=min(size(model.fullEV,2),numel(xValfullload0));
  y4=real(model.fullEV(k3,1:lastValue));
  yLabel='$\chi_1$';
  figname='_LAM36';
  dianame=strcat(modelfilename,figname);
  %set(0, 'DefaultFigureWindowState', 'minimized');
  figure(36);
  xPlot36=model.fulllambda(1:end)*xValfulllambda0Mult;
  [~,obergrenze]=markMins(xPlot36,y4);
  hold on
  if savefigures~=2
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   if k3==resEWs(1) && main.colorshift==0
    plot(model.fulllambda(1:end),zeros(size(model.fulllambda(1:end))),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
    grid on
   end
   
   plot(xPlot36,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);%,'Color',colJK ,'Color',colo);
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
   if model.savefigures==true
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')%dianame=strcat(modelfilename,figname);
   end
  end
  if savefigures==2
   cfig('wofuer')='TeX';
   cfig('yMax')=2.5e-4;%obergrenze;
   cfig('xMin')=lamMin;
   cfig('xMax')=lamMax;
   h8036=figure(8036);
   plotitJK(xPlot,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8036);
   title('(a)')
   markMins(xPlot,y4);
   print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))%dianame=strcat(modelfilename,figname);
  end
 end

 if ismember(37,plotfig)
  figure(37);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  end
  ls=res(1).stability_limit(end);
  x=xValload0;
  NrPlotVal=min(numel(x),NrLoad0);
  xlam=model.fulllambda(1:NrPlotVal);
  y4=transpose(real(model.fullEV(k3,1:NrPlotVal)))-(ls-xlam(:))/ls;
  plot(xlam,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
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

  %fignr=42
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
  Values=min(numel(xValload0)-1,numel(res(k3).NormR)-1);
  plot(xValload0(2:Values),(res(k3).NormR(2:Values)),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if savefigures==true
   dianame=strcat(modelfilename,'NormR_42');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 if ismember(43,plotfig) && numel(LAM)>0
  figure(43);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);

  hold on
  yData=res(k3).EVal(2:end);
  yNumel=numel(yData);
  xPlot43=lambda(2:yNumel+1)*xValfulllambda0Mult;
  plot(xPlot43,yData,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\chi$');
  grid on
  title(modelfilename)
  ylim([min(-1,max(real(LAM(end)),-0.2)),max(min(real(LAM(end)),2),1)])
  if savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_LAM43.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_LAM43.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_LAM43.pdf'),'-dpdf')
  end

 end
 
 if ismember(44,plotfig)
  figure(44);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  
  plot(lambda(3:end),(res(k3).EWd1l(3:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('$\frac{(EW(x+\Delta)-EW(x-\Delta)}{(2*\Delta)}$','Interpreter','latex');
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
  if savefigures==true
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
 
 %fignr=45
 if ismember(45,plotfig)
  figure(45);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  %x=model.fulllambda(1:end)*xValfulllambda0Mult; %xValfullload0
  y4=abs(model.fullEV(k3,1:numel(FullxPlot0)));
  plot(FullxPlot0,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   plot(FullxPlot0,0.*y4,'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  xlabel('$\lambda$','Interpreter','latex');
  yLabel='$|\chi|$';
  ylabel(yLabel,'Interpreter','latex');
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
   figname='_LAM45';
   dianame=strcat(modelfilename,figname);
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
  elseif savefigures==2
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('LineStyle')=MyLines{main.colorshift+1};
    cfig('XAxisLocation')='bottom';
    cfig('xMin')=lamMin;
    %cfig('xMax')=lamMax;
    cfig('yMax')=1;
    cfig('yMin')=0;
    h8015=figure(8015);
    plotitJK(xPlot0,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8015);
    if k3==resEWs(1) && main.colorshift==0
    elseif k3==resEWs(1) && main.colorshift==1
     legend('B32OS','B32OSH')
     dianame=strcat(modelfilename,figname);
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))%
    end
  else
   title(modelfilename)
  end
 end

 
 if ismember(46,plotfig)
  figure(46);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  x=model.fulllambda(1:end)*xValfulllambda0Mult; %xValfullload0
  y4=real(model.fullEV(k3,1:numel(x)))+imag(model.fullEV(k3,1:numel(x)));
  plot(x,y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   plot(x,zeros(size(x)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
   grid on
   grid minor
  end
  xlabel('$\lambda$','Interpreter','latex');
  yLabel='$Re(\chi)+imag(\chi)$';
  ylabel(yLabel,'Interpreter','latex');
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
   figname='_LAM46';
   dianame=strcat(modelfilename,figname);
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
  elseif savefigures==2
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('LineStyle')=MyLines{main.colorshift+1};
    cfig('XAxisLocation')='bottom';
    cfig('xMin')=lamMin;
    cfig('xMax')=lamMax;
    cfig('yMax')=1;
    cfig('yMin')=0;
    h8015=figure(8015);
    plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8015);
    if k3==resEWs(1) && main.colorshift==0
    elseif k3==resEWs(1) && main.colorshift==1
     legend('B32OS','B32OSH')
     dianame=strcat(modelfilename,figname);
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))%
    end
  else
   title(modelfilename)
  end
 end

 fignr=47;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$|A_0\cdot r_1|$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.NormeigvecA0r(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end


  fignr=48;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$|R1|$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.NormR1(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end


   fignr=49;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='r Kt0 r';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.rKt0rij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end

    fignr=50;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$r^{CT} Kt0 r$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.rCTKt0rij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end

    fignr=51;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$r^{NCT} \cdot Kt0 \cdot r$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.rNCTKt0rij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end

 fignr=52;
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=real(model.fullEV(k3,:));
  stepNRtmp=min(stepNR,numel(y4));
  zPlot=imag(model.fullEV(k3,1:stepNRtmp));
  if main.allowComplex==2
   y19imag=imag(model.fullEV(k3,1:numel(xPlot0)));
   y4(y19imag==0)=NaN;
   zPlot(y19imag==0)=NaN;
  elseif main.allowComplex==0
   y4(y19imag~=0)=NaN;
   zPlot(y19imag~=0)=NaN;
  end
  plot3(xPlot0(1:stepNRtmp),y4(1:stepNRtmp),zPlot,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
   daspect([1 1 1])
     view([1,1,1])
  title(modelfilename)
  if k3==resEWs(1) && main.colorshift==0
   camroll(120)
   grid on
   grid minor
  end
  yLabel='$Re(\chi_1)$';
  ylabel(yLabel,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  if strcmp(main.typeofanalysis,'KNL2')
   maxy4=max([y4,1.0]);
    bbb = gca();
    obergrenze=min(2,max(round(abs(maxy4),2,'significant'),maxy4));
    bbb.YLim = [-.5,obergrenze];
    cfig('yMax')=obergrenze;
  end
   figname=strcat('_LAM',num2str(fignr));
   dianame=strcat(modelfilename,figname);
  if savefigures==true
   set(gcf,'renderer','Painters')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  elseif savefigures==2
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('LineStyle')=MyLines{main.colorshift+1};
    cfig('XAxisLocation')='bottom';
    cfig('xMin')=lamMin;
    cfig('yMax')=1;
    cfig('yMin')=0;
    h8015=figure(8015);
    plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8015);
    if k3==resEWs(1) && main.colorshift==0
    elseif k3==resEWs(1) && main.colorshift==1
     legend('B32OS','B32OSH')
     dianame=strcat(modelfilename,figname);
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))%
    end
  else
   title(modelfilename)
  end
 end

     fignr=53;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$\cos(\varphi)=\frac{r_1\cdot A_0\cdot r_1}{|r_1|\,|A_0\,r_1|}$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.cosphirij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  %set(gca, 'YScale', 'log')%also enable the negative one
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  %plot(xPlot(1:Values),-yData(1:Values),'LineStyle','-','Marker','none','LineWidth',0.5,'Color',colJK);
  %xlim([-.65 -.6])
  %xlim([-.8 -.5])
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end

 fignr=54;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$\cos(\varphi)=\frac{r_1\cdot A_0\cdot r_1}{|r_1|\,|A_0\,r_1|}$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.cosphirij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  set(gca, 'YScale', 'log')%also enable the negative one
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  plot(xPlot(1:Values),-yData(1:Values),'LineStyle','-','Marker','none','LineWidth',0.5,'Color',colJK);
  %xlim([-.8 -.4])
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end
 
     fignr=55;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$Re(r) \cdot Kt0 \cdot Re(r)$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.RerNCTKt0Rerij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  %set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  dianame=strcat(modelfilename,'_normR',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
  end
 end
 
  
     fignr=56;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$Re(r) \cdot A_0 \cdot Re(r)$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.RerCTKt0Rerij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  %set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  bbb=gca;
  myylim=ylim;
  if myylim(2)>10
   bbb.YLim = [myylim(1),10];
   if myylim(1)<-10
    bbb.YLim = [-10,10];
   end
  end
  bbb.YLim = [-2*10^-7,4*10^-7];
  dianame=strcat(modelfilename,'_rAr',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
  end
 end
 
      fignr=57;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='$Re(r) \cdot A \cdot Re(r)$';
  ylabel(ylabelJK,'Interpreter','latex');
  if main.closall==true
   title(modelfilename,'Interpreter','none')
  end
  yData=model.RerARerij(:,k3);
  Values=min(numel(xPlot),numel(yData));%-3 because fulllamdba
  %set(gca, 'YScale', 'log')
  title(modelfilename,'Interpreter','none')
  plot(xPlot(1:Values),yData(1:Values),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  bbb=gca;
  myylim=ylim;
  if myylim(2)>10
   bbb.YLim = [myylim(1),10];
   if myylim(1)<-10
    bbb.YLim = [-10,10];
   end
  end
  bbb.YLim = [-10^-7,4*10^-7];
  dianame=strcat(modelfilename,'_rAr',num2str(fignr));
  if savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
  end
 end
 
 fignr=58;
  if ismember(fignr,plotfig)
   figure(fignr);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
   plot(lambda(2:end),res(k3).OC1(2:end),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
   xlabel('Lambda');
   ylabel('abs((r0*a) + (v*v))','Interpreter','latex');
   set(gca, 'YScale', 'log')
   if k3==resEWs(1) && main.colorshift==0
    grid on
    grid minor
   end
   if savefigures==true
    print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_X1.svg'))
    print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_X1.png'))
    print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'_X1.pdf'),'-fillpage')
   end
  end
 
 end %for k3
 
  if all(isnan(lambda))
   lambda=model.lambda0;
  end

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
  plot(lambda,zeros(size(lambda)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambda(1:lastvalue),relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  markMins(lambda,relDetKtx);
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\det(\mathbf K_T)$','Interpreter','latex');%/\det(\mathbf K_T)_0
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
  ylabel('$\sqrt[N]{\det(\mathbf K_T)}$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'sqrtDetKt901.pdf'),'-fillpage')
  end
 end
 
 idxs=[];
 if ismember(906,plotfig) || ismember(902,plotfig)
  h906=figure(906);
  yPlot=1./cell2mat(model.dxidl);
  idxsMin=islocalmin(yPlot,'FlatSelection', 'all');
  idxsMax=islocalmax(yPlot,'FlatSelection', 'all');
  idxs=idxsMin;
  idxs(idxsMax)=true;
  if ismember(906,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h906, 'Visible', 'off');
  end
  hold on
  xPlot906=model.lambda0*xValfulllambda0Mult;
  plot(xPlot906(1:numel(yPlot)),yPlot,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  xlabel('$\lambda$','Interpreter','latex');
  yLabel='$\mathrm d|P|/\mathrm du_{ave} [\textrm{N}/\textrm{m}]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  xPlot(idxs)
  grid on
  grid minor
  if main.savefigures>=true && ismember(906,plotfig)
   dianame=strcat(modelfilename,'dldi906');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    remove(cfig,{'yMax','yMin','xMin'});
    cfig('NaNs')='keep';
    cfig('xMin')=-1.5;
    cfig('xMax')=1.5;
    plotitJK(xPlotHalf(1:numel(yPlot)),yPlot,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,8906);
   end
  end
 end

 %plotfig=902;
 if ismember(902,plotfig)
  figure(902);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthuJK);
  %lambdaplot=model.load0;%lambda;
  lambdaplot=model.lambda0*xValfulllambda0Mult;%lambda;
  NrPlot=min(numel(y),find(~isnan(lambdaplot),1,'last'));
  %xPlot=lambdaplot(1:NrPlot);
  yPlot=y(1:NrPlot);
  if main.flipAxis
   plot(y(idxs),xPlot0(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
   plot(yPlot,xPlot0(1:NrPlot),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  else
   plot(xPlot0(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
   plot(xPlot0(1:NrPlot),yPlot,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  end
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
%   xLabel=xlabelload;%'$\lambda$';
  xlabel(myxlabelload,'Interpreter','latex');
  yLabel='$u_{ave} [\textrm{m}]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures>=true
   dianame=strcat(modelfilename,'xi902');
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    remove(cfig,{'yMax','yMin','xMin'});
    cfig('yMax')=0.08;
    cfig('NaNs')='keep';
    cfig('xMin')=lamMin;
    cfig('xMax')=max(lambdaplot);
    cfig('LineStyle')=MyLines{1};
    cfig('lineColor')=MyColours{1};
    h8902=figure(8902);
    if main.flipAxis
     plotitJK(yPlot,xPlot0(1:NrPlot),'Output/Figures/PlotIt/',yLabel,lambdaLabel,dianame,cfig,h8902);
     plot(y(idxs),lambdaplot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    else
     plotitJK(xPlot0(1:NrPlot),yPlot,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8902);
     plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    end
    print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
    print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
   end
  end
  %disp([xPlot,yPlot])
  %format longG
  %[~,lami]=InterpolateJK(1-xPlot,yPlot)%#ok<NOPRT,ASGLU> %fullEV=1-xPlot,lambda0=yPlot
 end

 if ismember(9021,plotfig)
  figure(9021);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthurHM);
  lambdaplot=model.load0;%lambda;
  NrPlot=min(numel(y),find(~isnan(lambdaplot),1,'last'));
  %plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  xPlot=lambdaplot(1:NrPlot);
  yPlot=y(1:NrPlot);
  plot(xPlot,yPlot,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='mean displacement $\xi_{common} [\textrm{m}]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
   grid on
   grid minor
  if main.savefigures>=true
   dianame=strcat(modelfilename,'xi902');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    cfig('yMax')=3;
    plotitJK(xPlot,yPlot,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,89021);
    %plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   end
  end
 end
 
 if ismember(9022,plotfig)
  figure(9022);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthuHM);
  lambdaplot=model.load0;%lambda;
  NrPlot=min(numel(y),find(~isnan(lambdaplot),1,'last'));
  %plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  xPlot=lambdaplot(1:NrPlot);
  yPlot=y(1:NrPlot);
  plot(xPlot,yPlot,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
%   xLabel=xlabelload;%'$\lambda$';
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='arc-length $\xi_{PH} [\textrm{m}]$';
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
    plotitJK(xPlot,yPlot,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,89022);
    %plot(lambdaplot(idxs),y(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
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
  xPlot=(lambdaplot(2:end)+lambdaplot(1:end-1))/2;
  lastValue=size(y,1);
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(xPlot(1:lastValue),y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
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
  lastValue=size(y4,1);
  lambdaplot=model.lambdainput;
  xPlot=(lambdaplot(2:lastValue+1)+lambdaplot(1:lastValue))/2;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(xPlot,y4(:,1),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
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
  y4=cell2mat(model.dxidl);
  lastValue=numel(y4);
  lambdaplot=xPlot0(1:lastValue);%model.lambdainput(1:lastValue);
  plot(lambdaplot,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d\xi(JK)/d\lambda [\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dxidl905.pdf'),'-fillpage')
  end
 end
 
 if ismember(907,plotfig)
  figure(907);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y=cell2mat(model.arclengthuJK);
  lambdaplot=model.lambda0;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
%   plot(lambdaplot,cell2mat(model.arclengthurHM),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
 % plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
%   plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(y,lambdaplot(1:numel(y)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  ylabel('$\lambda$','Interpreter','latex');
  xlabel('arc-length $\xi [\textrm{m}]$','Interpreter','latex');
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
  lambdaplot0=xPlot0;%model.load0;%=lambda;
  plot(lambdaplot0(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabel('$d\lambda/dv [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  idxs=markMins(lambdaplot0,y4);
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw909.pdf'),'-fillpage')
  end
 end
 
 if ismember(908,plotfig)
  figure(908);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=cell2mat(model.vMaxJK);
  lambdaplot=xPlot0;
  x=lambdaplot(1:numel(y4));
  if main.flipAxis
   plot(y4(idxs),x(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   plot(y4,x,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  else
   plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   plot(x,y4,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  end
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  %legend('arclengthurJK','arclengthuJK')
  %plot([0 0],[0 0],'LineWidth',eps(0));
  xlabel(myxlabelload,'Interpreter','latex');
  yLabel='$v_{max} [\textrm{m}]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  %if k3==resEWs(1) && main.colorshift==0
   grid on
   grid minor
  %end
  if main.savefigures>=true
   dianame=strcat(modelfilename,'w908');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   if main.savefigures==2
    h8908=figure(8908);
    cfig('xMin')=lamMin;
    cfig('LineStyle')=MyLines{1};
    cfig('lineColor')=MyColours{1};
    if main.flipAxis
     plotitJK(y4,x,'Output/Figures/PlotIt/',yLabel,lambdaLabel,dianame,cfig,h8908);
     plot(y4(idxs),x(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    else
     plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8908);
     plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    end
    %plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8908);
    title('(a)')
    print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
    print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
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
  plot(lambda,zeros(size(lambda)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(lambda(1:lastvalue),relDetKtx,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  plot([0 0],[0 0],'LineWidth',eps(0));
  plot(lambda(idxs),relDetKtx(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$\det(\mathbf K_T)$','Interpreter','latex');%/\det(\mathbf K_T)_0
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
  ylabel('$w_{max} [\textrm{m}]$','Interpreter','latex');
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
  %xPlot=xValload0;%lambda
  xPlot=model.fulllambda(1:end-3)*xValfulllambda0Mult;
  yPlot=relDetKtx;
  yLabel='$\det\mathbf K_T/\det(\mathbf K_T)_0$';
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
   cfig('xMin')=lamMin;
   cfig('xMax')=lamMax;
   h8913=figure(8913);
   hold on
   plotitJK(xPlot,yPlot,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8913);
   title('(b)')
   markMins(xPlot,yPlot);
   remove(cfig,{'yMax'});
   %markMins(xPlot,yPlot,gca());
   print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
  end
 end
 
 if ismember(914,plotfig)
  figure(914);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  xPlot=model.lambda0*xValfulllambda0Mult;
  lastvalue=min(numel(model.DetKtx),numel(xPlot));
  maxDet=max(model.DetKtx(1:lastvalue));
  relDetKtx=model.DetKtx(1:lastvalue)./maxDet;
  %xplot=lambda(1:lastvalue);
%   xplot=xValload0;
  y4=relDetKtx;
  plot(xPlot,zeros(size(xPlot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  plot(xPlot,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
%   xLabel='$\lambda$';
  xlabel(xlabelload,'Interpreter','latex');
  yLabel='$\displaystyle\frac{\det\mathbf  K_T}{{\scriptstyle \max}(\det\mathbf K_T)}$';
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
  %markMins(lambdaplot,dudl)
  if main.savefigures>=true
   figname='DetKt914';
   dianame=strcat(modelfilename,figname);
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   if main.savefigures==2
    remove(cfig,{'xMin'});
    cfig('NaNs')='keep';
    cfig('closing')=NaN;
    cfig('lineColor')=colJK;
    cfig('Faktor')=1;
    cfig('xMin')=-1.5;
    cfig('xMax')=1.5;
    cfig('yMin')=-.5;
    cfig('yMax')=1;
    cfig('XAxisLocation')='bottom';
    h8914=figure(8914);
    plotitJK(xPlot,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8914);
    %title('(b)')
    if k3==resEWs(1) && main.colorshift==0
     xline(0,'k','LineWidth',1)
     plot(xPlot,zeros(size(xPlot)),'Color',[0 0 0],'LineWidth',1)
     plotitJK(xPlot,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8914);
    end
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
  markMins(lambdaplot,y4);
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw915.pdf'),'-fillpage')
  end
 end

  if ismember(917,plotfig) || ismember(916,plotfig)
  h917=figure(917);
  if ismember(917,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h917, 'Visible', 'off');
  end
  hold on
  %y=cell2mat(model.arclengths);
  dudl=1./cell2mat(model.dudl);
  %lambdaplot=model.load0;%=lambda;
  lambdaplot=xPlot0;
  %plot(lambdaplot,zeros(size(lambdaplot)),'LineStyle','-','Marker','none','LineWidth',1,'Color','k')
  %plot(lambdaplot,y(:,5),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{1});
  %plot(lambdaplot,cell2mat(model.arclengthurJK),'LineStyle','--','Marker',markJK,'LineWidth',1.5,'Color',MyColours{2});
  %plot(lambdaplot,cell2mat(model.arclengthuHM),'LineStyle','-.','Marker',markJK,'LineWidth',1.5,'Color',MyColours{3});
  plot(lambdaplot(1:numel(dudl)),dudl,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabel('$d\lambda/dv [1/\textrm{m}]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  idxs917=markMins(lambdaplot,dudl);
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'dldw909.pdf'),'-fillpage')
  end
  %lambdaplot=model.load0;%=lambda
  %dudl=1./cell2mat(model.dudl);
  %figure(917);
  %idxs=markMins(lambdaplot,dudl);
 end
 
 if ismember(916,plotfig)
  figure(916);
  hold on
  y4=cell2mat(model.uMaxJK);
  x=xPlot0(1:numel(y4));
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  yLabel='$u_{max} [\textrm{m}]$';
  if main.flipAxis
   plot(y4(idxs917),x(idxs917),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   plot(y4,x,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
   xlabel(yLabel,'Interpreter','latex');
   ylabel(myxlabelload,'Interpreter','latex');
  else
   plot(x(idxs917),y4(idxs917),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   plot(x,y4,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
   xlabel(myxlabelload,'Interpreter','latex');
   ylabel(yLabel,'Interpreter','latex');
  end
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  title(modelfilename,'Interpreter','none')
   grid on
   grid minor
  if main.savefigures>=true
   dianame=strcat(modelfilename,'u916');
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   if main.savefigures==2
    cfig('xMin')=lamMin;
    cfig('xMax')=max(lambdaplot);
    cfig('LineStyle')=MyLines{1};
    cfig('lineColor')=MyColours{1};
    h8916=figure(8916);
    if main.flipAxis
     plotitJK(y4,x,'Output/Figures/PlotIt/',yLabel,lambdaLabel,dianame,cfig,h8916);
     plot(y4(idxs917),x(idxs917),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    else
     plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8916);
     plot(x(idxs917),y4(idxs917),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor',[1 1 1])
    end
    %plotitJK(x,y4,'Output/Figures/PlotIt/',lambdaLabel,yLabel,dianame,cfig,h8916);
    title('(b)')
    print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
    print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
   end
  end
 end

 if ismember(918,plotfig) ||  ismember(919,plotfig) || ismember(920,plotfig)|| ismember(969,plotfig)
  steps=numel(model.displacementsJK);
  [Dirs,Nodes]=size(model.displacementsJK{2});
  if Dirs<3
   error('MyPrgm:Dim','at least three directions')
   %assert(Dirs>=3,'at least three directions')
  end
  y1VecdisplacementsJK=NaN(steps,Nodes);
  y2Vec=NaN(steps,Nodes);
  y3Vec=NaN(steps,Nodes);
  tmp=max(6,floor(Nodes/2));
  NodesSchleife=2:tmp;
  for step=1:steps
   y1VecdisplacementsJK(step,:)=model.displacementsJK{step}(1,:);
   y2Vec(step,:)=model.displacementsJK{step}(2,:);
   y3Vec(step,:)=model.rotJK{step}(3,:);
  end
  if ismember(918,plotfig)
   figure(918);
   hold on
   x=model.load(1:steps);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   for node=NodesSchleife
    plot(x,y1VecdisplacementsJK(:,node),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',MyColours{mod(node,19)+1});
   end
   %plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(xlabelload,'Interpreter','latex');
   yLabel='$u [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'u918');
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   end
  end
  fignr=969;
  if ismember(fignr,plotfig)
   figure(fignr);
   hold on
   x=model.load(1:steps);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   for node=3
    plot(x,y1VecdisplacementsJK(:,node),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',MyColours{mod(node,19)+1});
    plot(x(idxs917),y1VecdisplacementsJK(idxs917,node),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   end
   %plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(xlabelload,'Interpreter','latex');
   yLabel='$u(l/2) [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'u',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   end
  end
  if ismember(919,plotfig)
   figure(919);
   hold on
   x=model.load(1:steps);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   for node=NodesSchleife
    plot(x,y2Vec(:,node),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',MyColours{mod(node,19)+1});
   end
   %plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(xlabelload,'Interpreter','latex');
   yLabel='$v [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'v919');
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   end
  end
  
  if ismember(920,plotfig)
   figure(920);
   hold on
   x=model.load(1:steps);
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   for node=NodesSchleife
    plot(x,y3Vec(:,node),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',MyColours{mod(node,19)+1});
   end
   %plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(xlabelload,'Interpreter','latex');
   yLabel='$w [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'w920');
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8920=figure(8920);
     %plotitJK(x,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,h8920);
     plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     plotitJK(x,y4,'Output/Figures/PlotIt/',xlabelload,yLabel,dianame,cfig,h8920);
     %print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'_.eps'))
    end
   end
  end
 end
 
  fignr=943;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\mathbf{r} \cdot \mathbf{\ddot{K}}_T \cdot \mathbf{r}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
   title(modelfilename,'Interpreter','none')
  bbb.YAxisLocation = 'origin';
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=model.rddotKtr(3:Xlength);
  plot(xValload(3:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

   fignr=944;% 
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='model.ZweirKtt';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.ZweirKtt));
  y4=model.ZweirKtt(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'ZweirKtt_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
   fignr=945;% (945)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='model.rdotKtr';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
   title(modelfilename,'Interpreter','none')
  bbb.YAxisLocation = 'origin';
  if strcmp(xBezug,'n')%lambda
   xPlot=lambda(1:end)*xValfulllambda0Mult; %xValfullload0
   xlabel('$\lambda$','Interpreter','latex');
  elseif strcmp(xBezug,'1')
   xPlot=lambda(1:end); %xValfullload0
   xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  else
   xPlot=xValload(1:end);%x=lambda(3:end)
   xlabel(xlabelload,'Interpreter','latex');
   xlim([1040000 1130000])
  end
  Xlength=min(numel(xValload),numel(model.rdotKtr));
  y4=model.rdotKtr(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 
   fignr=946;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='model.rKtr';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rKtr));
  y4=model.rKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
fignr=947;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\frac{ \ddot{\chi_1} }{ \mathbf{r}_1\cdot\mathbf{\ddot{K}}_T\cdot\mathbf{r}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=res(k3).EWd2l(1:Xlength)./model.rddotKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  myylim=ylim;
  if myylim(1)>0
   if myylim(2)>.9 && myylim(2)<1
    bbb.YLim = [0,1];
   else
    bbb.YLim = [0,myylim(2)];
   end
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'EWd2lprddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 fignr=948.1;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\frac{ -2 r_1 \dot{K}_T \dot{r}_1 }{ \mathbf{r}_1\cdot\mathbf{\ddot{K}}_T\cdot\mathbf{r}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=model.ZweirKtt(1:Xlength)./model.rddotKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'ZweirKttPrddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 fignr=948.5;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\frac{ -2 r_1 \dot{K}_T \dot{r}_1 }{ \mathbf{r}_1\cdot\mathbf{\ddot{K}}_T\cdot\mathbf{r}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=model.ZweirKtt(1:Xlength)./model.rddotKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'ZweirKttPrddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 fignr=949;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\frac{ 2 \dot{r}_1 [\dot{K}_T-\chi_1\,(\dot{K}_T)_0] \dot{r}_1 }{ \mathbf{r}_1\cdot\mathbf{\ddot{K}}_T\cdot\mathbf{r}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=2*model.t_KB1_t(1:Xlength)./model.rddotKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  myylim=ylim;
  if myylim(2)<1
   bbb.YLim = [myylim(1),1.01];
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'2t_KB1_tPrddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 fignr=950;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\frac{1}{\mathbf{r} \cdot \mathbf{\ddot{K}}_T \cdot \mathbf{r}}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=1./model.rddotKtr(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'OnerddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
    fignr=951;% 
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$r \dot{Kt} t$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.ZweirKtt));
  y4=(model.rdotKtt(1:Xlength));
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtt_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
  fignr=952;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='${\dot{r}_1 [\dot{K}_T-\chi_1\,(\dot{K}_T)_0] \dot{r}_1 }$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [-10,16];
  bbb.YLim = [-5e-5,5e-5];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=model.t_KB1_t(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'t_KB1_t_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
   fignr=953;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='${-2\cdot\dot{r}_1 [\dot{K}_T-\chi_1\,(\dot{K}_T)_0] \dot{r}_1 }$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
   title(modelfilename,'Interpreter','none')
  bbb.YAxisLocation = 'origin';
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=-2*model.t_KB1_t(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  ylim([-0.0001,0.0001])
  if model.savefigures==true
   dianame=strcat(modelfilename,'t_KB1_t_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
     fignr=954;% 
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$r \dot{Kt} t$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.ZweirKtt));
  y4=2*(model.rdotKtt(1:Xlength));
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtt_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
  fignr=955;% (43)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\displaystyle \frac{\mathbf{r}_1 \cdot \mathbf{\ddot{A}} \cdot \mathbf{r}_1}{\ddot{\chi}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=model.rddotKtr(1:Xlength)./res(k3).EWd2l(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  myylim=ylim;
  if myylim(1)>0
   if myylim(2)>.9 && myylim(2)<1
    bbb.YLim = [0,1];
   else
    bbb.YLim = [0,myylim(2)];
   end
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'rddotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
    fignr=956;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='$\displaystyle \frac{-2 \dot{r}_1 [\dot{A}-\chi_1\,(\dot{A})_0] \dot{r}_1 }{\ddot{\chi}_1}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xValload),numel(model.rddotKtr));
  y4=(-2*model.t_KB1_t(1:Xlength))./res(k3).EWd2l(1:Xlength);
  plot(xValload(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
   myylim=ylim;
  if myylim(2)<1
   bbb.YLim = [myylim(1),1.001];
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'t_KB1_t_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 if ismember(957,plotfig) && numel(lambda)>0
  fignr=957;
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);

  hold on
  lenMatch=numel(model.imagValues);
  xPlot=lambda(1:lenMatch-2)*xValfulllambda0Mult;
  plot(xPlot, model.imagValues(2:end-1),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  ylabel('number of imag EW');
  grid on
  title(modelfilename)
  if model.savefigures==true
   dianame=strcat(modelfilename,'_ImagNR',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

    fignr=958;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
   grid minor
  end
  if strcmp(xBezug,'n')%lambda
   xPlot=lambda(1:end)*xValfulllambda0Mult; %xValfullload0
   xlabel('$\lambda$','Interpreter','latex');
  elseif strcmp(xBezug,'1')
   xPlot=lambda(1:end); %xValfullload0
   xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  else
   xPlot=xValload(1:end);%x=lambda(3:end)
   xlabel(xlabelload,'Interpreter','latex');
   xlim([1040000 1130000])
  end
  ylabelJK='model.rKt0r';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xPlot),numel(model.rKt0r));
  y4=model.rKt0r(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
%   myylim=ylim;
%   if myylim(2)>10 && myylim(1)<10
%    bbb.YLim = [myylim(1),10];
%    if myylim(1)<-10
%     bbb.YLim = [-10,10];
%    end
%   end
  if model.savefigures==true
   dianame=strcat(modelfilename,'rKt0r_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

     fignr=959;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  if strcmp(xBezug,'n')%lambda
   xPlot=lambda(1:end)*xValfulllambda0Mult; %xValfullload0
   xlabel('$\lambda$','Interpreter','latex');
  elseif strcmp(xBezug,'1')
   xPlot=lambda(1:end); %xValfullload0
   xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  else
   xPlot=xValload(1:end);%x=lambda(3:end)
   xlabel(xlabelload,'Interpreter','latex');
   %xlim([1040000 1130000])
  end
  ylabelJK='model.NormKt0r';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xPlot),numel(model.NormKt0r));
  y4=model.NormKt0r(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'NormKt0r_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 fignr=960;
  if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=(model.d2ksi);
  lastValue=numel(y4);
  lambdaplot=xPlot0(1:lastValue);
  %xlim([.94 1.03])
  plot(lambdaplot,y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  xlabel('$\lambda$','Interpreter','latex');
  ylabel('$d^2\xi(JK) / (d\lambda)^2 [\textrm{m}^2]$','Interpreter','latex');
  title(modelfilename,'Interpreter','none')
  grid on
  grid minor
  if main.savefigures==true
   print('-dpdf',strcat('Output/Figures/PDF/',modelfilename,'d2xidl2_',num2str(fignr),'.pdf'),'-fillpage')
  end
  end

     fignr=961;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if k3==resEWs(1) && main.colorshift==0
   grid on
  end
  if strcmp(xBezug,'n')%lambda
   xPlot=lambda(1:end)*xValfulllambda0Mult; %xValfullload0
   xlabel('$\lambda$','Interpreter','latex');
  elseif strcmp(xBezug,'1')
   xPlot=lambda(1:end); %xValfullload0
   xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  else
   xPlot=xValload(1:end);%x=lambda(3:end)
   xlabel(xlabelload,'Interpreter','latex');
   xlim([1040000 1130000])
  end
  ylabelJK='model.RerKt0Imr';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xPlot),numel(model.RerKt0Imr));
  y4=model.RerKt0Imr(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'RerKt0Imr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end

 fignr=962;% ()
  if ismember(fignr,plotfig) || ismember(fignr+1,plotfig)
  h962=figure(fignr);
  if ismember(fignr,plotfig)
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  else
   set(h962, 'Visible', 'off');
  end
  hold on
  y4=1./cell2mat(model.dvdl);
  lambdaplot=model.load0;%=lambda;
  plot(lambdaplot(1:numel(y4)),y4,'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',MyColours{4});
  idxs=markMins(lambdaplot,y4);
 end
 
 fignr=963;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  y4=cell2mat(model.phiMaxJK);
  lambdaplot0=xPlot0;
  x=lambdaplot0(1:numel(y4));
  plot(x(idxs),y4(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
  plot(x,y4,'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5,'Color',colJK);
  if exist('xlimJK','var')
   xlim(xlimJK)
  end
  xlabel(myxlabelload,'Interpreter','latex');
  yLabel='$\phi_{max} [-]$';
  ylabel(yLabel,'Interpreter','latex');
  title(modelfilename,'Interpreter','none')
   grid on
   grid minor
  if main.savefigures>=true
   dianame=strcat(modelfilename,'phiMax',num2str(fignr));
   print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
  end
 end


 fignr=964;% ()
 if ismember(fignr,plotfig)
  steps=min(numel(model.uDiffJK),numel(xPlot0));
  [~,Nodes]=size(model.uDiffJK(2));
  y1VecuDiffJK=NaN(steps,Nodes);
  for step=2:steps
   y1VecuDiffJK(step,:)=model.uDiffJK(step);
  end
  yPlot=y1VecuDiffJK;
  [~,idxs]=max(model.duEdl);
   x=xPlot0(1:steps);
   figure(fignr);
   hold on
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   plot(x,yPlot(1:steps),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5);
   plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(myxlabelload,'Interpreter','latex');
   yLabel='$s [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'udiff',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8fignr=figure(8000+fignr);
     plotitJK(x,yPlot(1:steps),'Output/Figures/PlotIt/','$\lambda$',yLabel,dianame,cfig,h8fignr);
     plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     title('(b)')
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
     print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
    end
   end
 end


  fignr=965;% ()
 if ismember(fignr,plotfig)
  steps=min(numel(model.duEdl),numel(xPlot0));
  yPlot=model.duEdl;
  [~,idxs]=max(yPlot);
   x=xPlot0(1:steps);
   figure(fignr);
   hold on
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   plot(x,yPlot(1:steps),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5);
   xlabel(myxlabelload,'Interpreter','latex');
   yLabel='$\frac{Workincrement}{Loadincrement}=\frac{\partial W}{\partial P} [\frac{\textrm{N\,m}}{\textrm{N}}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'udiff',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8fignr=figure(8000+fignr);
     plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     plotitJK(x,yPlot(1:steps),'Output/Figures/PlotIt/','$\lambda$',yLabel,dianame,cfig,h8fignr);
     plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     title('(b)')
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
     print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
    end
   end
 end

   fignr=966;% ()
 if ismember(fignr,plotfig)
  steps=min(numel(model.duEdl),numel(xPlotHalf));
  %duEdl=([NaN;model.duEdl(1:end-1)]+model.duEdl)/2;
  lambdainputHalf=(model.lambdainput(1:end-1)+model.lambdainput(2:end))/2;
  yPlot=model.duEdl(1:end-1).*lambdainputHalf;
   x=xPlotHalf(1:steps);
  [~,xPoints,yHalfPoints,~] =markWende(x,yPlot);
   figure(fignr);
   hold on
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   plot(x,yPlot(1:steps),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5);
   plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(myxlabelload,'Interpreter','latex');
   yLabel='$\frac{\textrm{Workincrement}}{\textrm{Loadincrement}}=\frac{\mathrm dW}{\mathrm dP} [\frac{\textrm{N\,m}}{\textrm{N}}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'udiff',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8fignr=figure(8000+fignr);
     %plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     plotitJK(x,yPlot(1:steps),'Output/Figures/PlotIt/','$\lambda$',yLabel,dianame,cfig,h8fignr);
     plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     %title('(a)')
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
     print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
    end
   end
 end

    fignr=967;% ()
 if ismember(fignr,plotfig)
  steps=min(numel(model.duEdl(1:end-1)),numel(model.uDiffJK));
  %duEdl=([NaN;model.duEdl(1:end-1)]+model.duEdl)/2;
  lambdainputHalf=(model.lambdainput(1:end-1)+model.lambdainput(2:end))/2;
  yPlot=model.duEdl(1:end-1).*lambdainputHalf;
   x=model.uDiffJK(1:steps);
  [~,xPoints,yHalfPoints,~] =markWende(x,yPlot);
   figure(fignr);
   hold on
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   plot(x,yPlot(1:steps),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5);
   plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel('uDiffJK','Interpreter','latex');
   yLabel='$\frac{\textrm{Workincrement}}{\textrm{Loadincrement}}=\frac{\mathrm dW}{\mathrm dP} [\frac{\textrm{N\,m}}{\textrm{N}}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'udiff',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8fignr=figure(8000+fignr);
     %plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     plotitJK(x,yPlot(1:steps),'Output/Figures/PlotIt/','$s$',yLabel,dianame,cfig,h8fignr);
     plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     %title('(b)')
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
     print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
    end
   end
 end

    fignr=968;% ()
 if ismember(fignr,plotfig)
  steps=min(numel(model.duEdl),numel(xPlotHalf));
  %duEdl=([NaN;model.duEdl(1:end-1)]+model.duEdl)/2;
  %lambdainputHalf=(model.lambdainput(1:end-1)+model.lambdainput(2:end))/2;
  yPlot=model.duEdl(1:end-1);
   x=xPlotHalf(1:steps);
  %[~,xPoints,yHalfPoints,~] =markWende(x,yPlot);
   idxs=islocalmax(yPlot,'FlatSelection', 'all');
   figure(fignr);
   hold on
   set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
   plot(x,yPlot(1:steps),'LineStyle',lineStyleJK,'Marker',markJK,'LineWidth',1.5);
   plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
   xlabel(myxlabelload,'Interpreter','latex');
   yLabel='$\displaystyle\frac{\mathrm ds}{\mathrm d\lambda} [\textrm{m}]$';
   ylabel(yLabel,'Interpreter','latex');
   title(modelfilename,'Interpreter','none')
   grid on
   grid minor
   if main.savefigures>=true
    dianame=strcat(modelfilename,'udiff',num2str(fignr));
    print('-dpdf',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-fillpage')
    if main.savefigures==2
     h8fignr=figure(8000+fignr);
     %plot(xPoints,yHalfPoints,'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     plotitJK(x,yPlot(1:steps),'Output/Figures/PlotIt/','$\lambda$',yLabel,dianame,cfig,h8fignr);
     plot(x(idxs),yPlot(idxs),'o','Color',[0 0 0],'MarkerSize',12,'LineWidth',3,'MarkerFaceColor',[1 1 1])
     %title('(b)')
     print ('-depsc',strcat('Output/Figures/PlotIt/',dianame,'.eps'))
     print ('-dsvg',strcat('Output/Figures/PlotIt/',dianame,'.svg'))
    end
   end
 end
 
 fignr=970;% (970)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
   grid minor
  end
  xlabel(xlabelload,'Interpreter','latex');
  ylabelJK='model.RrA0Rr';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  title(modelfilename,'Interpreter','none')
  bbb.YAxisLocation = 'origin';
  Xlength=min(numel(xPlot),numel(model.rdotKtr));
  y4=model.RrA0Rr(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  myylim=ylim;
  if myylim(2)>10
   bbb.YLim = [myylim(1),10];
   if myylim(1)<-10
    bbb.YLim = [-10,10];
   end
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
 fignr=971;% (971)
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
   grid minor
  end
  xlabel(myxlabelload,'Interpreter','latex');
  ylabelJK='model.RrARr';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  title(modelfilename,'Interpreter','none')
  bbb.YAxisLocation = 'origin';
  Xlength=min(numel(xPlot),numel(model.rdotKtr));
  y4=model.RrARr(1:Xlength);
  plot(xPlot(1:Xlength),y4,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  myylim=ylim;
  if myylim(2)>10
   bbb.YLim = [myylim(1),10];
   if myylim(1)<-10
    bbb.YLim = [-10,10];
   end
  end
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtr_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
     fignr=972;% ()
 if ismember(fignr,plotfig)
  figure(fignr);
  set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',[FesterPosXNR(plotfig==get(gcf,'Number'))   FesterPosY   XBreite   YHohe]);
  hold on
  if main.colorshift==0
   grid on
   grid minor
  end
  yMult=1;
  if strcmp(xBezug,'n')%lambda
   xPlot=lambda(1:end)*xValfulllambda0Mult; %xValfullload0
   xlabel('$\lambda$','Interpreter','latex');
   yMult=1/xValfulllambda0Mult;
  elseif strcmp(xBezug,'1')
   xPlot=lambda(1:end); %xValfullload0
   xlabel('$\lambda gem Abaqus$','Interpreter','latex');
  else
   xPlot=xValload(1:end);%x=lambda(3:end)
   xlabel(xlabelload,'Interpreter','latex');
   xlim([1040000 1130000])
  end
  ylabelJK='$\frac{\mathbf{r}_\mathbf{u}\cdot\dot{\mathbf{K}}\cdot\mathbf{r}_\mathbf{u}}{\mathbf{r}_\mathbf{u}\cdot\mathbf{K}_0\cdot\mathbf{r}_\mathbf{u}} ; \frac{\hat{\mathbf{r}}^T\cdot\dot{\hat{\mathbf{K}}}\cdot\hat{\mathbf{r}}}{\hat{\mathbf{r}}^T\cdot\hat{\mathbf{K}}_0\cdot\hat{\mathbf{r}}}$';
  ylabel(ylabelJK,'Interpreter','latex');
  bbb = gca();
  %bbb.YLim = [0,1];
  %if main.closall==true
   title(modelfilename,'Interpreter','none')
  %end
  %bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  %xPlot=xValload;
  Xlength=min(numel(xPlot),numel(model.rKt0r));
  y4=model.rdotKtr(1:Xlength)./model.rKt0r(1:Xlength);
  plot(xPlot(1:Xlength),y4*yMult,'LineStyle','-','Marker','none','LineWidth',1.5,'Color',colJK);
  if model.savefigures==true
   dianame=strcat(modelfilename,'rdotKtr_p_rKt0r_',num2str(fignr));
   print('-dsvg',strcat('Output/Figures/SVG/',dianame,'.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',dianame,'.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',dianame,'.pdf'),'-dpdf')
  end
 end
 
end %function
