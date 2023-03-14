function [ fig ] = ...
 plotitJK(xVal,yMat,GFolder,xlabelJK,ylabelJK,diagramname,cfig,figh) 
 % plots the normalforce as a funktion of the cutingplace
 % Force...consists the force [ kcal / (mol*AA) ]
 % max... number of elements in force [ - ]
 % Folder... place wehre to put the graphs
 %set(0,'defaulttextinterpreter','latex')
 set(groot, 'DefaultAxesXColor', [0,0,0],'DefaultAxesYColor', [0,0,0],'DefaultAxesZColor', [0,0,0])%https://de.mathworks.com/matlabcentral/answers/465630-set-the-color-of-the-axes-per-default-to-black-instead-of-dark-gray#answer_377963

 if ~exist('GFolder','var')
  GFolder=pwd;
 end
 if ~strcmp(GFolder(end),'/')
  GFolder=strcat(GFolder,'/');
 end
 FigFolder=strcat(GFolder,"Figures/");
 Ort=pwd;
 cd(GFolder) %#ok<*MCCD>
 if isunix
 !mkdir -p Figures
 end
 cd(Ort)
 
 %#!
 %author:Johannes Kalliauer
 %university:TU Wien
 %created:2019-04-10

 %% == Description ==
 %

 %% == Last Changes ==
 %
  

 
 
 %% == Setting non existing input to default ==
 
 tmp=get(0,'defaultfigureposition');
 XBreite0=tmp(3);
 YHohe0=tmp(4);
 
   if ~exist('xlabelJK','var')
    XBeschriftung='x';
   else
    XBeschriftung=xlabelJK;
   end
   if ~exist('ylabelJK','var')
    ylabelJK='y';
   end
   if ~exist('diagramname','var')
    diagramname='Fig';cfig = containers.Map;
   end
   if ~exist('cfig','var')
    cfig = containers.Map;
   end
%    if ~exist('kuint','var')
%     kuint='';
%    end
%    if ~exist('order','var')
%     order=(1);
%    end
    if isKey(cfig,'order')
     order=cfig('order');
    else
     order=NaN;
    end
    if isKey(cfig,'kuint')
     kuint=cfig('kuint');
    else
     kuint='';
    end
%    if ~exist('xlog','var')
%     xlog=0;
%    end
   
   if ~exist('cfig','var')
    cfig = containers.Map;
   end
   
%    if exist('cfig','var')
    if isKey(cfig,'xlog')
     xlog=cfig('xlog');
    else
     xlog=false;
    end
    if isKey(cfig,'ylog')
     ylog=cfig('ylog');
    else
     ylog=false;
    end
    if isKey(cfig,'yMin')
     yMin=cfig('yMin');
     yMin0=yMin;
    else
     yMin=-Inf;
     yMin0=+1;%anything positiv
    end
    if isKey(cfig,'yMax')
     yMax=cfig('yMax');
     yMax0=yMax;
    else
     yMax=Inf;
     yMax0=-1;%anything negativ
     if verLessThan('matlab','9.7')
       warning('Mprgm:Matlab:too:Old','Maybe Matlab R2019a is too old')
     else
       if gca().YLim(2)~=1
         yMax=gca().YLim(2);
       end
     end
    end
    if isKey(cfig,'wofuer')
     wofuer=cfig('wofuer');
    else
     wofuer='plot';
    end
    if isKey(cfig,'xlabelmoveY')
     xlabelmoveY=cfig('xlabelmoveY');
    else
     xlabelmoveY=0;
    end
    if isKey(cfig,'Cxticks')
     Cxticks=cfig('Cxticks');
    end
    if isKey(cfig,'xticklabels')
     xticklabelsJK=cfig('xticklabels');
    end
    if isKey(cfig,'Cyticks')
     Cyticks=cfig('Cyticks');
    end
%     if isKey(cfig,'legend')
%      CLabels=cfig('legend');
%     end
    if isKey(cfig,'line') % replaced by LineStyle
     MyLineStyle='-';
     line=cfig('line');
     MyMarker='none';
     warning('MyPrgm:cfig','line should be replaced with cfig(LineStyle)')
    else
     MyLineStyle='-';
     MyMarker='o';
     %line='-o';
     %line=strcat(MyLineStyle,MyMarker);
    end
    if isKey(cfig,'MyMarker')
     MyMarker=cfig('MyMarker');
    end
    if isKey(cfig,'MarkerSize')
     MarkerSize=cfig('MarkerSize');
    else
     MarkerSize=get(0,'DefaultLineMarkerSize');
    end
    if isKey(cfig,'LineStyle')
     
     MyLineStyleMulti=cfig('LineStyle');
     if iscell(MyLineStyleMulti)
      MyLineStyle=MyLineStyleMulti{1};
      %lineMult=MyLineStyle;
      lineMultSet=true;
     else
      MyLineStyle=MyLineStyleMulti;
      lineMultSet=false;
     end
    end
    if isKey(cfig,'MyMarkerMulti')
     MyMarkerMulti=cfig('MyMarkerMulti');
     if ~isKey(cfig,'MyMarker')
      MyMarker=MyMarkerMulti{1};
     end
    end
     if isKey(cfig,'MyLineStyleMulti')
      MyLineStyleMulti=cfig('MyLineStyleMulti');
     end
    if isKey(cfig,'lineMult')
     %lineMult=cfig('lineMult');
     lineMultSet=true;
    else
     if ~exist('lineMultSet','var')
      lineMultSet=false;
     end
    end
    if isKey(cfig,'LineWidthMulti')
     LineWidthMulti=cfig('LineWidthMulti');
    end
    if isKey(cfig,'LineWidth')
     LineWidth=cfig('LineWidth');
     %LineWidthSet=true;
    else
     %LineWidthSet=false;
     LineWidth=2;
    end
    if isKey(cfig,'NaNs')
     NaNs=cfig('NaNs');
     if strcmp(NaNs,'keep') || strcmp(NaNs,'xkeep')
      xNaNs='keep';
     else
      xNaNs='del';
     end
     if strcmp(NaNs,'keep') || strcmp(NaNs,'ykeep')
      yNaNs='keep';
     else
      yNaNs='del';
     end
    else
     %NaNs='del';
     xNaNs='del';
     yNaNs='del';
    end
    if isKey(cfig,'XAxisLocation')
     XAxisLocation=cfig('XAxisLocation');
    else
     if isempty(yMat)
      return
     end
     yMaxWert=max(yMat(:));
     yMinWert=min(yMat(:));
     if max(yMaxWert,yMax0)>=0 && min(yMinWert,yMin0)<=0
      XAxisLocation='origin';
     else
      XAxisLocation='bottom';
     end
    end
    if isKey(cfig,'YAxisLocation')
     YAxisLocation=cfig('YAxisLocation');
    else
     YAxisLocation='left';
    end
    if isKey(cfig,'lineColor')
     MyColours=cfig('lineColor');
     if iscell(MyColours)
      MyColours1=MyColours{1};
     else
      MyColours1=MyColours;
     end
     if strcmp(MyColours,'JKextended')
      MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
      MyColours1=[0, 0.4470, 0.7410];
     end
     %lineColorSet=true;
    else
     MyColours={[0, 0.4470, 0.7410],	[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] 	,	[0.4660, 0.6740, 0.1880], 	[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],	[0, 0, 1],[0, 0.5, 0],	[1, 0, 0],	[0, 0.75, 0.75],[0.75, 0, 0.75],[0.75, 0.75, 0],[0.25, 0.25, 0.25],'y','m','c','g','k'};
     MyColours1=[0, 0.4470, 0.7410];
     %lineColorSet=false;
    end
    if isKey(cfig,'legend')
     MyLegendElements=cfig('legend');
    end
    if isKey(cfig,'xMin')
     xMin=cfig('xMin');
    else
     xMin=-Inf;
    end
    if isKey(cfig,'xMax')
     xMax=cfig('xMax');
    else
     xMax=Inf;
    end
    if isKey(cfig,'xLim')
     xLimJK=cfig('xLim');
     xMin=xLimJK(1);
     xMax=xLimJK(2);
    end
    if isKey(cfig,'yLim')
     yLimJK=cfig('yLim');
     yMin=yLimJK(1);
     yMax=yLimJK(2);
    end
    if isKey(cfig,'XBreite0')
     XBreite0=cfig('XBreite0');
    end
    if isKey(cfig,'YHohe0')
     YHohe0=cfig('YHohe0');
    end
    if isKey(cfig,'FontSize')
     FontSize0=cfig('FontSize');
    else
     FontSize0=(get(0,'defaultaxesfontsize')+get(0,'defaulttextfontsize'))/2;
    end
    if isKey(cfig,'fitmin')
     fitting=cfig('fitmin');
    else
     fitting=NaN;
    end
    if isKey(cfig,'closing')
     closing=cfig('closing');
    else
     closing='all';
    end
    if isKey(cfig,'xBorders')
     if strcmp(cfig('xBorders'),'exact')
      xMin=min(xVal(:));
      xMax=max(xVal(:));
     end
    end
    if isKey(cfig,'grid')
     gridJK=cfig('grid');
    else
     gridJK='on';%gridJK='minor';
    end
    if isKey(cfig,'Faktor')
     Faktor=cfig('Faktor');
    else
     Faktor=1.54;
    end
    if isKey(cfig,'outputformats')
     outputformats=cfig('outputformats');
    else
     outputformats='es'; %e..eps, s...svg
    end
    if isKey(cfig,'aspectratio')
      aspectratio=cfig('aspectratio');
%       if strcmp(aspectratio,'pbaspect([1 1 1])')
%         pbaspect([1 1 1])
%       elseif strcmp(aspectratio,'daspect')
%         daspect([1 1 1])
%       elseif strcmp(aspectratio,'equal')
%         axis equal
%       else
%         warning('MyPrgm:Unknown','aspectratio=%s unkown',aspectratio)
%       end
    end

    
   if strcmp(XAxisLocation,'origin')
    if all(yMat(:)<0) && yMax==Inf
     yMax=0;
    end
    if all(yMat(:)>0) && yMin==-Inf
     yMin=0;
    end
   end
  

   [Xzeilen,Xspalten]=size(xVal);
   if (Xzeilen>1 && Xspalten>1); warning('MyProgram:Dim','x is a matrix'); end
   if (Xzeilen==1 && Xspalten>1); xVal=transpose(xVal); [Xzeilen,Xspalten]=size(xVal); end
   assert(Xzeilen==1 || Xspalten==1,'x darf keine Matrix sein')
   [zeilen,spatlten]=size(yMat);
   if (zeilen~=Xzeilen && spatlten==Xzeilen); yMat=transpose(yMat);[zeilen,spatlten]=size(yMat); end
   if (zeilen==1 && spatlten>100 && Xzeilen>100 && Xspalten==1)
    error('MyProgram:Strange','Mabye transpose one input and change that inputvectors have the same length')
   end
   %yMat
   %size(yMat)
   if any(all(isnan(yMat))) %wenn es in min einer Spalte/Wertereihe nur NaN gibt
    %disp('Z138')
     if exist('MyLegendElements','var')
      LegendSize=size(MyLegendElements,1);
      toDelete=all(isnan(yMat), 1);
      DeleteElem=toDelete(1:LegendSize);
      MyLegendElements(DeleteElem) = [];
     end
    yMat(:,all(isnan(yMat), 1)) = [];
   end 
   %size(yMat)
   if size(yMat,2)==1
    NrYLines=1;
   else
    NrYLines=size(yMat,2);
   end
   
   if any(isnan(xVal)) && strcmp(xNaNs,'del')
    %disp('Z146')
    Werte=[xVal yMat];
    xVal(all(isnan(Werte), 2), :) = [];
    yMat(all(isnan(Werte), 2), :) = [];
   end
   %disp(any(isnan(yMat)));
   if all(any(isnan(yMat))) && strcmp(yNaNs,'del') %wenn es in jeder Spalte immer min einen NaN gibt
    %disp('Z138')
    xVal(all(isnan(yMat), 2)) = [];
    yMat(all(isnan(yMat), 2), :) = [];

   end
   
   %% ==Calculation ==
   if NrYLines==1 && ~isnan(order)
    xValMittel=xVal;%-median(xVal);
    if any(yMat)
     PoR2 = polyfitn(double(xValMittel),double(yMat),double(order));
     %PoR2.Coefficients
    else
     PoR2.R2=NaN;
    end
    [q2m,S] = polyfit(double(xValMittel),double(yMat),double(order));
    %q2m
    xValUnique=unique(xVal);
    [y_fit,delta] = polyval(q2m,xValUnique,S); %#ok<ASGLU>
   end
if all(kuint~=0) && ~isnan(order) && NrYLines==1
 SpNBeschrm=strcat('k =~',int2str(q2m(1)*single(factorial(order))),kuint,' (R2=',num2str(PoR2.R2),')');
 disp(SpNBeschrm);
else
 SpNBeschrm='';
end
 

 cd(FigFolder)
 if strcmp(closing,'all')%default value
  close all
 end
 %% Figuredef
FesterPosY =uint16(300);
FesterPosX =uint16(558);
FontName = 'times';
%Faktor=1.54;
XBreite=(XBreite0)*Faktor; %416
YHohe=(YHohe0)*Faktor; %312
FontSize=FontSize0*Faktor;

gcfPostition=[FesterPosX   FesterPosY   XBreite   YHohe]; 
if strcmp(wofuer,'TeX')
 tmp1=.2;
 tmp2=.15;
 if FontSize>=17.5
  tmp2=.16;
 end
 tmp3=.98-tmp1;
 %tmp4=.99-tmp2;
 tmp4=.94-tmp2;
 gcaPosition=[tmp1 tmp2 tmp3 tmp4];
 Zusatzinfos=false;
elseif strcmp(wofuer,'PNG')
 gcaPosition=[0.20 0.15 0.79 0.83];
 %Position=[0.01 0.15 0.79 0.83];
 Zusatzinfos=false;
else
 gcaPosition=[0.10 0.11 0.88 0.89]; %[0.10 0.08 0.88 0.90]; %416x312 %[0.08 0.08 0.92 0.91]
 Zusatzinfos=true;
end
 
if usejava('jvm')
 if exist('figh','var') && ~isempty(figh)
  figure(figh)
 else
  figure;
 end
 set(gca,'FontSize',FontSize*Faktor,'FontName',FontName,'Position',gcaPosition,'XAxisLocation',XAxisLocation,'YAxisLocation',YAxisLocation) % %[0.20 0.20 0.67 0.675])
 if xlog==1
  %semilogx(xVal,yMat)
  set(gca, 'XScale', 'log')
 end
 if ylog==1
  %   semilogy(xVal,yMat)
  set(gca, 'YScale', 'log')
 end
 set(gcf,'PaperUnits','points','PaperPositionMode','auto','PaperOrientation','landscape','Position',gcfPostition); %,'Position',gcfPostition
 box on
 hold on
 if exist('aspectratio','var')
      if strcmp(aspectratio,'pbaspect([1 1 1])')
        pbaspect([1 1 1])
      elseif strcmp(aspectratio,'daspect')
        daspect([1 1 1])
      elseif strcmp(aspectratio,'equal')
        axis equal
      else
        warning('MyPrgm:Unknown','aspectratio=%s unkown',aspectratio)
      end
 end
 if strcmp(gridJK,'off')
  grid off
 elseif strcmp(gridJK,'on')
  grid on
 elseif  strcmp(gridJK,'minor')
  grid on
  grid minor
 else
  warning('MyProgram:Unknown','grid not defined')
  grid on
 end
 
%set(0, 'DefaultFigureWindowState', 'minimized');
   if NrYLines==1

    if ~isnan(all(MyMarker)) && ~isnan(all(MyColours1))
     if ~strcmp(MyMarker,'none')
      plot(xVal,yMat,'Marker',MyMarker,'Color',MyColours1,'LineStyle',MyLineStyle,'MarkerSize',MarkerSize,'LineWidth',LineWidth)
     else
      plot(xVal,yMat,'Marker',MyMarker,'Color',MyColours1,'LineStyle',MyLineStyle,'LineWidth',LineWidth)
     end
    else
     plot(xVal,yMat,line,'LineStyle',MyLineStyle)
    end
    if isreal(xMin) && ~isinf(xMin) && ~isinf(xMax) && ~isnan(xMin) && ~isnan(xMax)
      if xMax>xMin
        xlim([xMin xMax])
      end
    end
    if NrYLines==1 && ~isnan(order)
     r2m=refcurve(q2m);
     %r2m.Color=[0 1 0];
     r2m.Color=MyColours1;
     %r2mn=refcurve(PoR2.Coefficients);r2mn.Color=[1 0 0];
     if order==2 && ~all(isnan(fitting))
      Diffs=yMat-(q2m(1)*xVal.^2+q2m(2)*xVal+q2m(3));
      Versatz=abs(min(Diffs));
      q2mU=q2m+[0,0,Versatz];
      q2mL=q2m-[0,0,Versatz];
      r2mU=refcurve(q2mU);r2mU.Color='r';
      r2mL=refcurve(q2mL);r2mL.Color='r';
      Auswahl=(Diffs<Versatz);
      if ~all(Auswahl)
       del=xVal(~Auswahl) %#ok<NOPRT,NASGU>
       plotitJK(xVal(Auswahl),yMat(Auswahl),GFolder,xlabelJK,ylabelJK,diagramname,kuint,order,xlog,cfig) 
      end
     end
     %plot(xValUnique,y_fit,'r-')
     %plot(xValUnique,y_fit+delta,'m-.',y_fit-delta,'m-.')
     %plot(xValUnique,y_fit+2*delta,'m--',xValUnique,y_fit-2*delta,'m--')
    end
    if isreal(yMin) && ~isnan(yMin) && isfinite(yMin)
     %disp(yMin)
     if yMax~=Inf && yMin<yMax
      ylim([yMin yMax])
     end
    end
   else %NrYLines>0
    for i=1:NrYLines
     xVali=xVal;
     yMati=yMat(:,i);
     %Werte=[xVal yMati];
     if any(isnan(yMati)) && strcmp(yNaNs,'del')
      %disp('Z211')
      xVali(any(isnan(yMati), 2), :) = [];
      yMati(any(isnan(yMati), 2), :) = [];
     end
     if any(isnan(xVal)) && strcmp(xNaNs,'del')
      disp('220')
      xVali(any(isnan(xVal), 2), :) = [];
      yMati(any(isnan(xVal), 2), :) = [];
     end
%      if lineColorSet==true
      %i
      %xVali
      %yMati
      if lineMultSet==true
       if exist('MyMarkerMulti','var')
        MyMarker=cell2mat(MyMarkerMulti(i));
       else
        MyMarker='none';
       end
       MyLineStyle=cell2mat(MyLineStyleMulti(i));
       if exist('LineWidthMulti','Var')
        LineWidth=cell2mat(LineWidthMulti(i));
       end%      elseif lineMultSet==true && numel(lineMult)>1
%       cell2mat(lineMult(i))
%       plot(xVali,yMati,cell2mat(lineMult(i)))
%      else
%       plot(xVali,yMati,'LineStyle',MyLineStyle)
%      end

      end
      plot(xVali,yMati,'LineStyle',MyLineStyle,'Marker',MyMarker,'Color',cell2mat(MyColours(i)),'LineWidth',LineWidth,'MarkerSize',MarkerSize)
%      elseif lineMultSet==true && numel(lineMult)>1
%       cell2mat(lineMult(i))
%       plot(xVali,yMati,cell2mat(lineMult(i)))
%      else
%       plot(xVali,yMati,'LineStyle',MyLineStyle)
%      end
     if isreal(yMin) && ~isnan(yMin) && isfinite(yMin)
      ylim([yMin yMax]) ;
     end
     if isreal(xMin) && ~isnan(xMin)
      assert(xMin<xMax,'xMin must be smaler than xMax')
      xlim([xMin xMax])
     end 
     
    end
    if exist('cfig','var')
     if exist('MyLegendElements','var')
      l=legend(MyLegendElements);   set(l,'Interpreter','latex')
      if isKey(cfig,'legendLocation')
       set(l, 'Location',cfig('legendLocation'))   %set(l, 'Position', [.32 0.75 0 0])
      end
     end
    end
   end
      if isKey(cfig,'xline')
       xline(cfig('xline'),'Color',[0 0 0],'LineWidth',2);
      end
      if isKey(cfig,'yline')
       %yl=yline(cfig('yline'),'Color',[0 0 0],'LineWidth',2);
       %yl.Color = [0 0 0];
       plot(xlim,[0 0],'Color',[0 0 0],'LineWidth',2)
      end

   xlh=xlabel(XBeschriftung,'Interpreter','latex');
   xlh.Position(2) = xlh.Position(2) + xlabelmoveY;
   if exist('Cxticks','var')
    xticks(Cxticks)
    if exist('xticklabelsJK','var')
     set(gca,'TickLabelInterpreter', 'latex');
     xticklabels(xticklabelsJK)
    end
   end
   if exist('Cyticks','var')
    yticks(Cyticks)
   end
   %ylh=ylabel(ylabelJK,'Interpreter','latex','Rotation',0); %#ok<NASGU>
   ylh=ylabel(ylabelJK,'Interpreter','latex'); %#ok<NASGU>
   %ylh.Position(1)=ylh.Position(1)-0.33;
   %t=title('Force');set(t,'Interpreter','latex')
   
   %ax.XAxisLocation = 'origin';
   
   if Zusatzinfos
    text(0.1,.95,{strcat(GFolder,"Figures/",diagramname),date},'Units','normalized','Interpreter','none'); %#ok<*DATE> % %#ok<DATE> 
    text(.5,.75,SpNBeschrm,'Color','blue','Units','normalized','Interpreter','none');
   end
   
   fig = gcf;

 %hold off
 %print ('-depsc','mini.eps')
 %!epstopdf "mini.eps" --outfile=mini.pdf
 %!./JK.sh
 if contains(outputformats,'e')
   print ('-depsc',strcat(diagramname,'.eps'))
 end

 %print('-fillpage',diagramname,'-dpdf')
 %orient portrait
 if contains(outputformats,'p')
   print('-dpng',strcat(diagramname,'.png'))
 end
 if contains(outputformats,'s')
   print('-dsvg',strcat(diagramname,'.svg'))
 end
 if contains(outputformats,'w') && ispc
  print('-dmeta',strcat(diagramname,'.emf'))%Enhanced metafile (Windows only)
 end

end

 
%dlmwrite(strcat(FigFolder,'Data.csv'), [xVal;yMat], 'delimiter', ',', 'precision', 7); 
%dlmwrite(strcat(FigFolder,'DataT.csv'), [xVal;yMat]', 'delimiter', ',', 'precision', 7); 

cd(Ort)

%close all force;


end % function
