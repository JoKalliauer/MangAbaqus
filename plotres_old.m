function plotres(res,model,plotfig,colJK,markJK)
%% Plot Graphs
%author:Michal Malendowski (2019-2020) and Johannes Kalliauer(2020-2021),
%later replaced by plotresMulti.m
 % 2...rho
 % 1...chi
 % 3..._velocity
 % 4_tanacceleration
 % 5_noracceleration
 % 6_totacceleration
 % 7_torque
 % 8_SinCos
 % 8.1 Sin/Cos
 % 11_Eigenwert ylim([0,1])
 % % 11.1 ylim detail
 % 12_real(Eigenwert) kein ylim
 % 13_abs(eigenwert)
 % 14_rho2
 
 if ~exist('model','class')
  modelfilename=model.filename;
 else
  modelfilename='unknown';
 end
 if ~exist('plotfig','var')
  plotfig=[1,2];
 end
 %MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};


 
 
 %lambda1 = res.lambdaorg;
 lambda = res.lambda;
 %V = res.V;
 %A = res.A;
 S = res.S;
 A0 = res.A0;%...total accerlation
 At = res.At;
 An = res.An;
 RHO = res.RHO;
 %RHO2 = res.RHO2;
 TAU = res.TAU;
 cosmu = res.cosmu;
 sinpsi = res.sinpsi;
 %       OC1 = OC1(Orth);
 %OCeig = res.OCeig;
 LAM = res.LAM;
 %R = res.R;
 %POS = res.POS;
 
 %OC0 = res.OC0;
 %OC1 = res.OC1;
 %OC2 = res.OC2;
 %OC3 = res.OC3;
 
 %coplanar = res.coplanar;
 
 %colo = [rand(), rand(), rand()];
 %colo = [0 0 0];
 
 kl = res.stability_limit(1); %#ok<NASGU>
 %stability_limit = res.stability_limit(2);
 
 
 %maxlam = max(lambda);
 maxlam = 1.0;
 lambda = lambda/maxlam;
 
 if model.savefigures==true
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
 
 if ismember(2,plotfig)
  figure(2);
  hold on
  plot(lambda(2:end),RHO(2:end),'LineStyle','-','Marker','o','LineWidth',1.5,'Color',colJK);%,'Color',colo);
  
  %        plot(lambda(2:end),RHO2(2:end),'LineStyle','none','Marker','.','Color',colo);
  %        plot(lambda,RHO3,'LineStyle','none','Marker','o','Color',colo);
  %        plot(lambda,RHO4,'LineStyle','--','Marker','x','Color',colo);
  xlabel('lambda');
  ylabel('rho');
  title('rho');
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
 
 lambda = reshape(lambda,length(lambda),1); LAM = reshape(LAM,length(lambda),1);
 LAM = LAM/maxlam;
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
  grid on
  %        legend('chi','eigenvalue','lambda','zero')
  title(modelfilename)
  xticks(0:10)
  yticks(-40:50)
  %ylim([0 min(20,lambda(end)+real(LAM(end)))])
  daspect([1 1 1])
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_chi.svg'))
  print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_chi.png'))
  print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_chi.pdf'),'-dpdf')
  end
 end

 
 if ismember(3,plotfig)
  figure(3)
  hold on
  plot(lambda(2:end),S(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('velocity');
  title('velocity');
  grid on
  nopeaks = (S(S<10*mean(S(:))));
  bbb = gca();
  bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_velocity.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_velocity.png'))
  end
 end
 
 if ismember(4,plotfig)
  figure(4)
  hold on
  At = real(At);
  plot(lambda(2:end),At(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('tangential acceleration');
  title('tangential acceleration');
  grid on
  nopeaks = (At(abs(At)<10*mean(abs(At(:)))));
  bbb = gca();
  bbb.YLim = [min([nopeaks;-eps(1)]),max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_tanacceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_tanacceleration.png'))
  end
 end
 
 if ismember(5,plotfig)
  figure(5)
  hold on
  plot(lambda(2:end),An(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('normal acceleration');
  title('normal acceleration');
  grid on
  nopeaks = (An(An<10*mean(An(:))));
  bbb = gca();
  bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_noracceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_noracceleration.png'))
  end
 end
 
 if ismember(6,plotfig)
  figure(6)
  hold on
  plot(lambda(2:end),A0(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  xlabel('lambda');
  ylabel('total acceleration');
  title('total acceleration');
  grid on
  nopeaks = (A0(A0<10*mean(A0(:))));
  bbb = gca();
  bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_totacceleration.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_totacceleration.png'))
  end
 end
 
 if ismember(7,plotfig)
  figure(7)
  hold on
  plot(lambda(2:end),abs(TAU(2:end)),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  %       plot(lambda,abs(TAU2),'LineStyle','--','Marker','o','Color',colo);
  xlabel('lambda');
  ylabel('torque');
  title('torque');
  grid on
  nopeaks = (TAU(TAU<10*mean(TAU(:))));
  bbb = gca();
  bbb.YLim = [0.0,max([nopeaks;eps(1)])];
  if model.savefigures==true
  print('-dsvg',strcat('Output/Figures/',modelfilename,'_torque.svg'))
  print('-dpng',strcat('Output/Figures/',modelfilename,'_torque.png'))
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
 
 % % %        figure(9)
 % % %        hold on
 % % %        plot(lambda(2:end),asin(sinpsi(2:end))*180/pi,'LineStyle','-','LineWidth',1.5,'Marker','x','Color',col);
 % % %        plot(lambda(2:end),acos(cosmu(2:end))*180/pi,'LineStyle','--','LineWidth',1.5,'Marker','o','Color',col);
 % % % %        plot(lambda(2:end),(asin(sinpsi(2:end))+acos(cosmu(2:end)))*180/pi,'LineStyle','-','Marker','.','Color',[1,0,0]);
 % % % %        plot(lambda(2:end),abs(asin(sinpsi(2:end))-acos(cosmu(2:end)))*180/pi,'LineStyle','--','Marker','.','Color',[1,0,0]);
 % % %        xlabel('Lambda');
 % % %        ylabel('angle');
 % % % %        legend('psi','nu','psi+nu','abs(psi-nu)')
 % % %        grid on
 % % % %        nopeaks = (TAU(TAU<10*mean(TAU(:))));
 % % %        bbb = gca();
 % % %        maxy = max(bbb.YTick);
 % % %        bbb.YTick = 0.0:45:maxy;
 
 % % %        figure(10)
 % % %        if (asin(sinpsi(2))*180/pi>=0)&&(acos(cosmu(2))*180/pi<=90)
 % % %            psi1 = asin(sinpsi)*180/pi;
 % % %            nu1 = acos(cosmu)*180/pi;
 % % %            psi2 = asin(sinpsi)*180/pi;
 % % %            nu2 = 360 - acos(cosmu)*180/pi;
 % % %            psi3 = 180 - asin(sinpsi)*180/pi;
 % % %            nu3 = acos(cosmu)*180/pi;
 % % %            psi4 = 180 - asin(sinpsi)*180/pi;
 % % %            nu4 = 360 - acos(cosmu)*180/pi;
 % % %        elseif (asin(sinpsi(2))*180/pi>=0)&&(acos(cosmu(2))*180/pi>=90)
 % % %            psi1 = asin(sinpsi)*180/pi;
 % % %            nu1 = acos(cosmu)*180/pi;
 % % %            psi2 = asin(sinpsi)*180/pi;
 % % %            nu2 = 360 - acos(cosmu)*180/pi;
 % % %            psi3 = 180 - asin(sinpsi)*180/pi;
 % % %            nu3 = acos(cosmu)*180/pi;
 % % %            psi4 = 180 - asin(sinpsi)*180/pi;
 % % %            nu4 = 360 - acos(cosmu)*180/pi;
 % % %        elseif (asin(sinpsi(2))*180/pi<=0)&&(acos(cosmu(2))*180/pi<=90)
 % % %            psi1 = 180 - asin(sinpsi)*180/pi;
 % % %            nu1 = acos(cosmu)*180/pi;
 % % %            psi2 = 180 - asin(sinpsi)*180/pi;
 % % %            nu2 = 360 - acos(cosmu)*180/pi;
 % % %            psi3 = 360 + asin(sinpsi)*180/pi;
 % % %            nu3 = acos(cosmu)*180/pi;
 % % %            psi4 = 360 + asin(sinpsi)*180/pi;
 % % %            nu4 = 360 - acos(cosmu)*180/pi;
 % % %        elseif (asin(sinpsi(2))*180/pi<=0)&&(acos(cosmu(2))*180/pi>=90)
 % % %            psi1 = 180 - asin(sinpsi)*180/pi;
 % % %            nu1 = acos(cosmu)*180/pi;
 % % %            psi2 = 180 - asin(sinpsi)*180/pi;
 % % %            nu2 = 360 - acos(cosmu)*180/pi;
 % % %            psi3 = 360 + asin(sinpsi)*180/pi;
 % % %            nu3 = acos(cosmu)*180/pi;
 % % %            psi4 = 360 + asin(sinpsi)*180/pi;
 % % %            nu4 = 360 - acos(cosmu)*180/pi;
 % % %        end
 % % %        subplot(2,2,1)
 % % %        hold on
 % % %        plot(lambda(2:end),psi1(2:end),'LineStyle','-','Marker','x','Color',colo);
 % % %        plot(lambda(2:end),nu1(2:end),'LineStyle','--','Marker','o','Color',colo);
 % % %        plot(lambda(2:end),psi1(2:end)+nu1(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
 % % %        bbb = gca();
 % % %        bbb.YLim = [0,360];
 % % %        bbb.YTick = 0.0:90:360;
 % % %        grid on
 % % %
 % % %        subplot(2,2,2)
 % % %        hold on
 % % %        plot(lambda(2:end),psi2(2:end),'LineStyle','-','Marker','x','Color',colo);
 % % %        plot(lambda(2:end),nu2(2:end),'LineStyle','--','Marker','o','Color',colo);
 % % %        plot(lambda(2:end),psi2(2:end)+nu2(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
 % % %        bbb = gca();
 % % %        bbb.YLim = [0,360];
 % % %        bbb.YTick = 0.0:90:360;
 % % %        grid on
 % % %
 % % %        subplot(2,2,3)
 % % %        hold on
 % % %        plot(lambda(2:end),psi3(2:end),'LineStyle','-','Marker','x','Color',colo);
 % % %        plot(lambda(2:end),nu3(2:end),'LineStyle','--','Marker','o','Color',colo);
 % % %        plot(lambda(2:end),psi3(2:end)+nu3(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
 % % %        bbb = gca();
 % % %        bbb.YLim = [0,360];
 % % %        bbb.YTick = 0.0:90:360;
 % % %        grid on
 % % %
 % % %        subplot(2,2,4)
 % % %        hold on
 % % %        plot(lambda(2:end),psi4(2:end),'LineStyle','-','Marker','x','Color',colo);
 % % %        plot(lambda(2:end),nu4(2:end),'LineStyle','--','Marker','o','Color',colo);
 % % %        plot(lambda(2:end),psi4(2:end)+nu4(2:end),'LineStyle','-','Marker','.','Color',[1,0,0]);
 % % %        bbb = gca();
 % % %        bbb.YLim = [0,360];
 % % %        bbb.YTick = 0.0:90:360;
 % % %        grid on
 % % %
 % % %        figure(3000)
 % % %        hold on
 % % %        plot(lambda1,abs(OC0),'LineStyle','--','Marker','s','Color',colo);
 % % %        hold on
 % % %        plot(lambda1,abs(OC1),'LineStyle','--','Marker','o','Color',colo);
 % % %        plot(lambda1,abs(OC2),'LineStyle','--','Marker','.','Color',colo);
 % % %        plot(lambda1,abs(OC3),'LineStyle','--','Marker','^','Color',colo);
 % % %        plot(lambda,1-OCeig,'LineStyle','-','Marker','x','Color',colo);
 % % %        plotaxes = gca();
 % % %        plotaxes.YScale = 'log';
 % % %        xlabel('Lambda');
 % % %        ylabel('Orthogonality checks');
 % % %     %    legend({'n*t','1 - (r_i * r_{i+1})'})
 % % %        legend({'n*t','T*N','T*B','N*B','1 - (r_i * r_s)'})
 % % %        grid on
 % % %        hold off
 % % %
 % % %        figure(3001)
 % % %        hold on
 % % %        plot(lambda(2:end),OCeig(2:end),'LineStyle','-','Marker','none','LineWidth',1.5,'Color',col);
 % % %        plot(lambda(2:end),OC0(2:end),'LineStyle','--','Marker','none','LineWidth',1.5,'Color',col);
 % % %        plotaxes = gca();
 % % % %        plotaxes.YScale = 'log';
 % % %        xlabel('Lambda');
 % % %     %    legend({'n*t','1 - (r_i * r_{i+1})'})
 % % % %        legend({'r * r'''})
 % % %        grid on
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
  disp([lambda(2:end), LAM(2:end)])
 end
 
 if ismember(12,plotfig)
  figure(12);
  hold on
  %        plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none');%,'Color',colo);
  %        plot(lambda(2:end), LAM(2:end),'LineStyle','--','Marker','none');%,'Color',colo);
  aaa = gca();
  if size(aaa.Children,1)<2
   %plot(lambda(1:end),lambda(1:end),'Color',[0 0 0],'LineStyle','--');
   %plot(lambda(1:end),0*lambda(1:end),'Color',[0 0 0]);
  end
  %plot(lambda(2:end), lambda(2:end)+LAM(2:end),'LineStyle','-','Marker','none','LineWidth',1.5);%,'Color',col,'Color',colo);
  plot(lambda(2:end),real(LAM(2:end)),'LineStyle','-','Marker',markJK,'LineWidth',1.5,'Color',colJK);%,'Color',colo);
  xlabel('Lambda');
  %        ylabel('Lambda + chi | or | chi');
  ylabel('real(EW)=real(chi-lambda)');
  grid on
  %yticks(-50:2:10)
  %        legend('chi','eigenvalue','lambda','zero')
  %ylim([-Inf max(min(20,real(LAM(end))),1)])
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
  hold on
  grid on
  grid minor
  xlabel('lambda');
  ylabel('rho2');
  bbb = gca();
  bbb.YLim = [0.0,1];
  bbb.XLim = [0 inf];
  title(modelfilename)
  bbb.XAxisLocation = 'origin';
  bbb.YAxisLocation = 'origin';
  plot(lambda(2:end),res.RHO2(2:end),'LineStyle','-','Marker','o','LineWidth',1.5,'Color',colJK);
  %yticks(0:.05:1)
  grid(bbb,'minor','on')
  if model.savefigures==true
   print('-dsvg',strcat('Output/Figures/SVG/',modelfilename,'_rho.svg'))
   print('-dpng',strcat('Output/Figures/PNG/',modelfilename,'_rho.png'))
   print('-fillpage',strcat('Output/Figures/PDF/',modelfilename,'_rho.pdf'),'-dpdf')
  end
  grid minor
 end
 if ismember(-14,plotfig)
  disp([lambda(2:end),res.RHO2(2:end)])
 end
 
end
