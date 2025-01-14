function [filename,lambda,BC,Nodes,Elements,load,dofpNode,rpLeft,leftnodes,rpRight,rightnodes] ...
= Kreis_2024(numofelm,lambda,RefLast,eltype,AbaqusRunsFolder,modelprops)
%university:TU Wien
%author: Johannes Kalliauer(2022-2023)
%uses code from TL_arch3D by Michal Malendowski (2019-2020), Johannes Kalliauer(2020-2023)
%created:2022

%% Input
% ~ ... length (not used, since it is received via modelprops)
% numofelm ... number of elements
% lambda ... load steps
% loadFactor .. Load factor to multiply the load
% eltype ... type of the Element
% AbaqusRunsFolder .. where is the folder to run Abaqus
% modelprops ... input-properties such as the youngs modulus or the length

%% Output 
% filename ... name of the abaqus-file
% lambda ... load steps (same as input?)
% BC ... boundary conditions
% Nodes .. nodenumber and coordinates
% Elements ... which nodes to which element
% load .. load-values for each lambda
% dofpNode ... number of dof per Node
% rpLeft .. left support (currently not used)
% leftnodes ... left support (currently not used)
% rpRight ... right support (currently not used)
% rightnodes ... right support (currently not used)

%% Last changes
% 2023-03-13 JK: removed useless razem -Output, added comments, added modelprops.length-check

%% Default input and input-check
if nargin<1
    lambda = 0.01:0.01:1;
%     lambda = 0.1;
end
if lambda(1) == 0
    lambda(1) = [];
end
if nargin<2
    numofelm = 50;
end
if length(eltype)>=5
 if strcmp(eltype(1:5),'B32OS')
  warning('MyProgram:Input','TL_arch3D only with rect tested')
 end
end

assert(~isempty(modelprops.length),'modelprops.length is []');

elmtype=eltype;

 MV=modelprops.MeterValue;



 
%% Span Vorgaben:
 R=modelprops.length*MV;% [m] 19.074
 W=215/360;% [-] Umdrehungen
 if sum(strcmp(fieldnames(modelprops.profil), 'E')) || sum(strcmp(fieldnames(modelprops.profil), 'nu'))
  Emodul=modelprops.profil.E;
 else
  warning('MyPrgm:Outdated','assuming E=2e+11N/m²');
  Emodul=2e+11/MV; % [N/m^2]
  nu=.25;
 end
 if sum(strcmp(fieldnames(modelprops.profil), 'h')) || sum(strcmp(fieldnames(modelprops.profil), 'b'))
  h=modelprops.profil.h;
  b=modelprops.profil.b;
 else
  warning('MyPrgm:Outdated','assuming h=0.2m and b=0.1m');
  h = 20e-2*MV;
  b = 10e-2*MV;
 end

 

 filename = ['Kreis_2024-',eltype,'-',num2str(numofelm(end)),'-f',num2str(RefLast),'-eps',num2str(modelprops.epsilon),'-u',num2str(MV),'-h',num2str(h),'-b',num2str(b),'-R',num2str(R)];
 
 %EA=Emodul*b*h;
 %EI=Emodul*b*h^3/12;
 %kappa=6/5;
 %GATilde=EA/(2*(1+nu)*kappa)

%% Load
 %p = -83.3*10^3*10^2*loadFactor; % [N/m ?]ca
 p=0; % [N/m ?]ca
 PJK=RefLast;%[N]
 load=-lambda*PJK;
  
%% Imperfection size
 %impsize = 0.001;
 impsize = 0;

%% Finite Elements Size
 resolution = W/double(numofelm(end));
 
%% Finite Element Model
 elmPerLength = round(W/(2*resolution))*2;

 alcoords = linspace(-W*pi,W*pi,elmPerLength+1)';
 alcoords(abs(alcoords)<1e-12) = 0;
 xcoords = R*sin(alcoords);
 ycoords = R*cos(alcoords);

%   plot(xcoords,ycoords,'mo-'); hold off
    
 Nodes = [ctranspose(1:length(alcoords)), xcoords, ycoords];
 Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
 
 rpLeft = Nodes(1,1);
 leftnodes = Nodes(1,1);
 TopNode=Nodes((elmPerLength/2)+1,1);
 rpRight = Nodes(end,1);
 rightnodes = Nodes(end,1);
 
 if strcmpi('B32',eltype(1:3))
  alcoords2=0.5*(alcoords(1:end-1,1:end) + alcoords(2:end,1:end));
  xcoords2 = R*sin(alcoords2);
  ycoords2 = R*cos(alcoords2);
  %zcoords2 = 0*alcoords2;
  Nodes2 = [Nodes(end,1) + ctranspose(1:size(alcoords2,1)),xcoords2,ycoords2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 end
  Nodes = [Nodes, zeros(size(Nodes,1),1)];
  
 PMal = zeros(size(Nodes,1),1);
 Pd = zeros(size(Elements,1),1);
 for i = 1:size(Elements,1)
     node1 = Elements(i,2); node2 = Elements(i,3);
     if strcmpi(elmtype(1:3),'B32')
         node3 = Elements(i,4);
         coords = [Nodes(node1,2:3); Nodes(node2,2:3); Nodes(node3,2:3)];
            % extract position on XY plane
           x1 = coords(1,1); x2 = coords(2,1); x3 = coords(3,1); 
           y1 = coords(1,2); %y2 = coords(2,2); 
           y3 = coords(3,2);

           Lload1 = abs(x2 - x1);
           Lload2 = abs(x3 - x2);

           PMal(node1) = PMal(node1) + 0.5*p*Lload1;
           PMal(node2) = PMal(node2) + 0.5*p*Lload1;
           PMal(node2) = PMal(node2) + 0.5*p*Lload2;
           PMal(node3) = PMal(node3) + 0.5*p*Lload2;
           
           Pd(i) = p*abs(x3 - x1)/sqrt((x3-x1)^2 + (y3-y1)^2);
     else
         coords = [Nodes(node1,2:3); Nodes(node2,2:3)];
            % extract position on XY plane
           x1 = coords(1,1); x2 = coords(2,1); 
           y1 = coords(1,2); y2 = coords(2,2); 

           Lload1 = abs(x2 - x1);

           PMal(node1) = PMal(node1) + 0.5*p*Lload1;
           PMal(node2) = PMal(node2) + 0.5*p*Lload1;
           
           Pd(i) = p*abs(x2 - x1)/sqrt((x2-x1)^2 + (y2-y1)^2);
     end
 end
 
 %P = [Nodes(:,1),P];
 %Pd = [Elements(:,1),Pd];
 
%  plotMesh(Nodes,Elements);
%  hold off

if ~exist(AbaqusRunsFolder, 'dir')
 if isunix
     mkdir(AbaqusRunsFolder);
 end
 if ispc
  warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
  return
 end
end
     
u1 = fopen([AbaqusRunsFolder,filename,'.inpData'],'w');

fprintf(u1,'*Part, name=Part-1\n');
fprintf(u1,'*Node\n');
fprintf(u1,'%d, %f, %f, %f \n',Nodes');
fprintf(u1,['*Element, type=',elmtype,'\n']);
if strcmpi(elmtype(1:3),'B32')
    fprintf(u1,'%d, %d, %d, %d\n',Elements');
else
    fprintf(u1,'%d, %d, %d\n',Elements');
end
fprintf(u1,'*Elset, elset=AllElements\n');
fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(:,1));
if length(Elements(:,1))/16~=floor(length(Elements(:,1))/16)
    fprintf(u1,'\n');
end

fprintf(u1,'** Section: Section-1  Profile: Profile-1\n');
fprintf(u1,'*Beam Section, elset=AllElements, material=Material-1, temperature=VALUES, section=RECT\n');
profile = [b, h];
fprintf(u1,'%f, %f \n',profile);
fprintf(u1,'0.,0.,-1\n');
fprintf(u1,'*End Part\n');

fprintf(u1,'*Assembly, name=Assembly\n');
fprintf(u1,'*Instance, name=Part-1-1, part=Part-1\n');
fprintf(u1,'*End Instance\n');

fprintf(u1,'*Nset, nset=allnodes, instance=Part-1-1, generate\n');
fprintf(u1,['1,',num2str(Nodes(end,1)),',1\n']);
fprintf(u1,'*Elset, elset=allelements, instance=Part-1-1, generate\n');
fprintf(u1,['1,',num2str(Elements(end,1)),',1\n']);

fprintf(u1,'*Nset, nset=leftend, instance=Part-1-1\n');
fprintf(u1,[num2str(rpLeft),'\n']);
fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
fprintf(u1,[num2str(rpRight),'\n']);
fprintf(u1,'*Nset, nset=TopNode, instance=Part-1-1\n');
fprintf(u1,[num2str(TopNode),'\n']);

fprintf(u1,'*End Assembly\n');

fprintf(u1,'*Material, name=Material-1\n');
fprintf(u1,'*Elastic\n');
fprintf(u1,[num2str(Emodul),', ',num2str(nu),'\n']);

%% Boundary conditions
  if strcmp(elmtype,'B32OSH') || strcmp(elmtype,'B31OSH') %%B31H,B33H 
   dofpNode=7;
  elseif  strcmp(elmtype,'B32') || strcmp(elmtype,'B31') || strcmp(elmtype,'B33') || strcmp(elmtype,'B32H') || strcmp(elmtype,'B31H') || strcmp(elmtype,'B33H')
   dofpNode=6;
  else
%    dofpNode=6;
   dofpNode=[];
  end
 BC = [dofpNode*(rpLeft - 1) + 1, 0
       dofpNode*(rpLeft - 1) + 2, 0;
       dofpNode*(rpRight - 1) + 1, 0;
       dofpNode*(rpRight - 1) + 2, 0;
       dofpNode*(Nodes(:,1)-1)+3, zeros(size(Nodes,1),1)
       dofpNode*(Nodes(:,1)-1)+4, zeros(size(Nodes,1),1)
       dofpNode*(Nodes(:,1)-1)+5, zeros(size(Nodes,1),1)];
  
  BCrpLeft=1:6;
  for i=1:numel(BCrpLeft)
   BC=[BC;dofpNode*(rpLeft - 1) + BCrpLeft(i), 0]; %#ok<AGROW>
  end
  [~,idx]=unique(BC(:,1));
  BC=BC(idx,:);
%%

fprintf(u1,'*Boundary\n');
fprintf(u1,'leftend,  1, 1\n');
fprintf(u1,'leftend,  2, 2\n');
fprintf(u1,'leftend,  1, 6\n');
fprintf(u1,'rightend,  1, 1\n');
fprintf(u1,'rightend,  2, 2\n');
fprintf(u1,'allnodes,  3, 3\n');
fprintf(u1,'allnodes,  4, 4\n');
fprintf(u1,'allnodes,  5, 5\n');


%impsize=0;

u3 = fopen([AbaqusRunsFolder,filename,'.inp'],'w');

if impsize~=0
    fprintf(u3,['*Imperfection, File=',filename,'-imperfections, STEP=1\n']);
    fprintf(u3,['1, ',num2str(impsize),'\n']);
end

fprintf(u3,['*Include, input=',filename,'.inpData\n']);
fprintf(u3,'** ----------------------------------------------------------------\n');
fprintf(u3,'*STEP, name=Lambda-1\n');
fprintf(u3,'*MATRIX GENERATE, STIFFNESS\n');
fprintf(u3,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
fprintf(u3,'*END STEP\n');
% % % 
% % % fprintf(u3,'** ----------------------------------------------------------------\n');
% % % fprintf(u3,'*Step, name=Eig-1, nlgeom=NO, perturbation\n');
% % % fprintf(u3,'*Buckle, eigensolver=lanczos\n');
% % % fprintf(u3,'10, , , , \n');
% % % 
% % %     fprintf(u3,'*Boundary, op=NEW, load case=1\n');
% % %     fprintf(u3,'leftend,  1, 1\n');
% % %     fprintf(u3,'leftend,  2, 2\n');
% % %     fprintf(u3,'rightend,  1, 1\n');
% % %     fprintf(u3,'rightend,  2, 2\n');
% % %     fprintf(u3,'*Boundary, op=NEW, load case=2\n');
% % %     fprintf(u3,'leftend,  1, 1\n');
% % %     fprintf(u3,'leftend,  2, 2\n');
% % %     fprintf(u3,'rightend,  1, 1\n');
% % %     fprintf(u3,'rightend,  2, 2\n');
% % %        
% % %     fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
% % %     fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), Pd(:,2)]');
% % %     
% % %     fprintf(u3,'*NODE PRINT, nset=allnodes\n');
% % %     fprintf(u3,'U\n');
% % %     
% % % fprintf(u3,'*END STEP\n');

stepnum = 1;
for k = 1:length(lambda)
    stepnum = stepnum + 1;
    fprintf(u3,'** ----------------------------------------------------------------\n');
    fprintf(u3,'** \n');
    fprintf(u3,['** STEP: Step-',num2str(stepnum),'\n']);
    fprintf(u3,'** \n');
    fprintf(u3,['*Step, name=Step-',num2str(stepnum),', nlgeom=YES\n']);
    fprintf(u3,'*Static\n');
    if k==1
       fprintf(u3,[num2str(lambda(k)),', ',num2str(lambda(k)),', ',num2str(lambda(k)*0.0001),', ',num2str(lambda(k)),'\n']);
    else
       fprintf(u3,[num2str(lambda(k)-lambda(k-1)),', ',num2str(lambda(k)-lambda(k-1)),', ',num2str((lambda(k)-lambda(k-1))*0.0001),', ',num2str(lambda(k)-lambda(k-1)),'\n']);
    end
    fprintf(u3,'** \n');
    fprintf(u3,'** LOADS\n');
    fprintf(u3,'** \n');

%     if k==1
%         fprintf(u3,'*Cload, OP=NEW\n');
%     else
%         fprintf(u3,'*Cload, OP=MOD\n');
%     end    
%     fprintf(u2,'Part-1-1.%d, 2, %f\n',[P(:,1), lambda(k)*P(:,2)]');

    if k==1
        fprintf(u3,'*Cload, OP=NEW\n');
    else
        fprintf(u3,'*Cload, OP=MOD\n');
    end    
    fprintf(u3,'TopNode, 2, %f\n',-lambda(k)*PJK);
    
%     if k==1
%          fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
%     else
%          fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=MOD\n');
%     end
%     fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), lambda(k)*Pd(:,2)]');
        
    fprintf(u3,'** \n');
    fprintf(u3,'** OUTPUT REQUESTS\n');
    fprintf(u3,'** \n');
    fprintf(u3,'*Restart, write, frequency=0\n');
    fprintf(u3,'** \n');
    fprintf(u3,'** FIELD OUTPUT: F-Output-1\n');
    fprintf(u3,'** \n');
    %fprintf(u3,'*Output, field, variable=PRESELECT\n');
    fprintf(u3,'*Output, field\n');
    fprintf(u3,'*Node Output\n');
    fprintf(u3,'CF, RF, TF, U\n');
    fprintf(u3,'*Element Output, directions=YES\n');
    fprintf(u3,'ESF1, LE, NFORC, NFORCSO, PE, PEEQ, PEMAG, S, SF\n');
    fprintf(u3,'** \n');
    fprintf(u3,'** HISTORY OUTPUT: H-Output-1\n');
    fprintf(u3,'** \n');
    fprintf(u3,'*Output, history, variable=PRESELECT\n');
    fprintf(u3,'*NODE PRINT, nset=allnodes\n');
    fprintf(u3,'U\n');
    fprintf(u3,'*EL PRINT\n');
    fprintf(u3,'\n');
    fprintf(u3,'SE\n');
    fprintf(u3,'*EL PRINT\n');
    fprintf(u3,'\n');
    fprintf(u3,'SF\n');
    fprintf(u3,'*End Step\n');
    
    stepnum = stepnum + 1;
    fprintf(u3,'** ----------------------------------------------------------------\n');
    fprintf(u3,['*STEP, name=Lambda-',num2str(stepnum),'\n']);
    fprintf(u3,'*MATRIX GENERATE, STIFFNESS\n');
    fprintf(u3,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
    fprintf(u3,'*END STEP\n');

end

fclose(u3);


end

