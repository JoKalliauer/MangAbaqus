function [filename,lambda,BC,Nodes,Elements,rpLeft,leftnodes,rpRight,rightnodes,razem] = TL_arch3D_Hinge(dummy,numofelm,lambda,loadFactor,eltype)
razem = [];
if nargin<1
    lambda = 0.05:0.05:0.35;
    numofelm = 50;
%     lambda = 0.1;
end
if lambda(1) == 0
    lambda(1) = [];
end

% pure SI units: Newtons, meters, Pascals, etc.
filename = ['TL_arch3D_Hinge-',eltype,'-',num2str(numofelm(end)),'-loadfac-',num2str(loadFactor)];

%% RECT
 h = 20e-2;
 b = 10e-2;
  
%% Span
 L = 600e-2; % m
 H = 240e-2; % m

%% Load
 p = -83.3*10^3*10^2*loadFactor;
  
%% Imperfection size
 impsize = 0.001;
 impsize = 0;

%% Finite Elements Size
 resolution = L/2/(numofelm(end));
 
%% Finite Element Model
 elmPerLength = round(L/resolution);

%%
    function f = arclength(x,L,H)
       f = -1/4*(L-2*x)*sqrt((16*H^2*(L-2*x)^2+L^4)/L^4) - (L^2*asinh((4*H*(L-2*x))/L^2))/16/H;
       x = 0;
       f = f - (-1/4*(L-2*x)*sqrt((16*H^2*(L-2*x)^2+L^4)/L^4) - (L^2*asinh((4*H*(L-2*x))/L^2))/16/H);
    end
 al_tot = arclength(L,L,H);
 DL = al_tot/elmPerLength;
 
 xcoords = zeros(elmPerLength/2+1,1);
  dx = zeros(elmPerLength/2,1);
 for i = 2:(elmPerLength/2+1)
    dx0 =  DL/sqrt(1 + 16*H^2/L^4*(L - 2*xcoords(i-1))^2);
    dx(i-1) = dx0;
    xcoords(i) = xcoords(i-1) + dx0;
 end
%  xcoords = linspace(0,L,elmPerLength+1)';
 xcoords = xcoords*((L/2)/max(xcoords));
 xcoords = [xcoords; L - xcoords(end-1:-1:1)];
 
 xcoords = 1e-12*round(1e12*xcoords);
 xcoords(abs(xcoords)<1e-12) = 0;
 ycoords = 4*H/L^2*xcoords.*(L - xcoords);
%   plot(xcoords,ycoords,'mo-'); hold off
    
 Nodes = [[1:length(xcoords)]', xcoords, ycoords];
 Elements = [[1:size(Nodes,1)-1]',Nodes(1:end-1,1),Nodes(2:end,1)];
 
 rpLeft = Nodes(1,1);
 leftnodes = Nodes(1,1);
 rpRight = Nodes(end,1);
 rightnodes = Nodes(end,1);
 
 if strcmpi(elmtype(1:3),'B32')
  xcoords2 = 0.5*(Nodes(1:end-1,2) + Nodes(2:end,2));
  ycoords2 = 4*H/L^2*xcoords2.*(L - xcoords2);
     
  Nodes2 = [xcoords2, ycoords2];
  Nodes2 = [Nodes(end,1) + [1:size(Nodes2,1)]',Nodes2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 end
 
  midnode1 = Nodes(Nodes(:,2)==L/2,1);
  Nodes = [Nodes; Nodes(midnode1,:)];
  midnode2 = size(Nodes,1);
  Nodes(end,1) = midnode2;
  Elements(Elements(:,2)==midnode1,2) = midnode2;
  Nodes = [Nodes, zeros(size(Nodes,1),1)];
 
 P = zeros(size(Nodes,1),1);
 Pd = zeros(size(Elements,1),1);
 for i = 1:size(Elements,1)
     node1 = Elements(i,2); node2 = Elements(i,3);
     if strcmpi(elmtype(1:3),'B32')
         node3 = Elements(i,4);
         coords = [Nodes(node1,2:3); Nodes(node2,2:3); Nodes(node3,2:3)];
            % extract position on XY plane
           x1 = coords(1,1); x2 = coords(2,1); x3 = coords(3,1); 
           y1 = coords(1,2); y2 = coords(2,2); y3 = coords(3,2);

           Lload1 = abs(x2 - x1);
           Lload2 = abs(x3 - x2);

           P(node1) = P(node1) + 0.5*p*Lload1;
           P(node2) = P(node2) + 0.5*p*Lload1;
           P(node2) = P(node2) + 0.5*p*Lload2;
           P(node3) = P(node3) + 0.5*p*Lload2;
           
           Pd(i) = p*abs(x3 - x1)/sqrt((x3-x1)^2 + (y3-y1)^2);
     else
         coords = [Nodes(node1,2:3); Nodes(node2,2:3)];
            % extract position on XY plane
           x1 = coords(1,1); x2 = coords(2,1); 
           y1 = coords(1,2); y2 = coords(2,2); 

           Lload1 = abs(x2 - x1);

           P(node1) = P(node1) + 0.5*p*Lload1;
           P(node2) = P(node2) + 0.5*p*Lload1;
           
           Pd(i) = p*abs(x2 - x1)/sqrt((x2-x1)^2 + (y2-y1)^2);
     end
 end
 
 P = [Nodes(:,1),P];
 Pd = [Elements(:,1),Pd];
 
%  plotMesh(Nodes,Elements);
%  hold off

u1 = fopen(['AbaqusRuns/',filename,'-model.inp'],'w');

fprintf(u1,'*Part, name=Part-1\n');
fprintf(u1,'*Node\n');
fprintf(u1,'%d, %f, %f, %f \n',Nodes');
fprintf(u1,['*Element, type=',elmtype,'\n']);
if strcmpi(elmtype(1:3),'B32')
    fprintf(u1,'%d, %d, %d, %d\n',Elements');
else
    fprintf(u1,'%d, %d, %d\n',Elements');
end
% fprintf(u1,'*Element, type=B32\n');
% fprintf(u1,'%d, %d, %d, %d\n',Elements');
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

fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+1,midnode1,midnode2]);
fprintf(u1,'*Connector Section, Elset=connector\n');
fprintf(u1,'JOIN\n');

fprintf(u1,'*End Assembly\n');

fprintf(u1,'*Material, name=Material-1\n');
fprintf(u1,'*Elastic\n');
fprintf(u1,'2e+11, 0.3\n');

%% Boundary conditions
 BC = [6*(rpLeft - 1) + 1, 0
       6*(rpLeft - 1) + 2, 0;
       6*(rpRight - 1) + 1, 0;
       6*(rpRight - 1) + 2, 0;
       6*(Nodes(:,1)-1)+3, zeros(size(Nodes,1),1)
       6*(Nodes(:,1)-1)+4, zeros(size(Nodes,1),1)
       6*(Nodes(:,1)-1)+5, zeros(size(Nodes,1),1)];
%%

fprintf(u1,'*Boundary\n');
fprintf(u1,'leftend,  1, 1\n');
fprintf(u1,'leftend,  2, 2\n');
fprintf(u1,'rightend,  1, 1\n');
fprintf(u1,'rightend,  2, 2\n');
fprintf(u1,'allnodes,  3, 3\n');
fprintf(u1,'allnodes,  4, 4\n');
fprintf(u1,'allnodes,  5, 5\n');

u2 = fopen(['AbaqusRuns/',filename,'-imperfections.inp'],'w');
fprintf(u2,['*Include, input=',filename,'-model.inp\n']);
fprintf(u2,'** ----------------------------------------------------------------\n');
fprintf(u2,'** \n');
fprintf(u2,'** STEP: imperfections\n');
fprintf(u2,'** \n');
fprintf(u2,'*Step, name=imperfections, nlgeom=NO, perturbation\n');
fprintf(u2,'*Buckle, eigensolver=lanczos\n');
fprintf(u2,'10, 0., , , \n');
fprintf(u2,'** \n');
fprintf(u2,'** LOADS\n');
fprintf(u2,'** \n');

% fprintf(u2,'*Cload, OP=NEW\n');
% fprintf(u2,'Part-1-1.%d, 2, %f\n',P');
fprintf(u2,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
fprintf(u2,'Part-1-1.%d, PY, %f\n',Pd');

fprintf(u2,'** \n');
fprintf(u2,'** OUTPUT REQUESTS\n');
fprintf(u2,'** \n');
fprintf(u2,'*Restart, write, frequency=0\n');
fprintf(u2,'** \n');
fprintf(u2,'** FIELD OUTPUT: F-Output-2\n');
fprintf(u2,'** \n');
fprintf(u2,'*Output, field, variable=PRESELECT\n');
fprintf(u2,'*NODE FILE\n');
fprintf(u2,'U\n');
fprintf(u2,'*End Step\n');

fclose(u2);

u3 = fopen(['AbaqusRuns/',filename,'.inp'],'w');

if impsize~=0
    fprintf(u3,['*Imperfection, File=',filename,'-imperfections, STEP=1\n']);
    fprintf(u3,['1, ',num2str(impsize),'\n']);
end

fprintf(u3,['*Include, input=',filename,'-model.inp\n']);
fprintf(u3,'** ----------------------------------------------------------------\n');
fprintf(u3,'*STEP, name=Lambda-1\n');
fprintf(u3,'*MATRIX GENERATE, STIFFNESS\n');
fprintf(u3,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
fprintf(u3,'*END STEP\n');

% fprintf(u3,'** ----------------------------------------------------------------\n');
% fprintf(u3,'*Step, name=Eig-1, nlgeom=NO, perturbation\n');
% fprintf(u3,'*Buckle, eigensolver=lanczos\n');
% fprintf(u3,'10, , , , \n');
% 
%        
%     fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
%     fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), Pd(:,2)]');
%     
%     fprintf(u3,'*NODE PRINT, nset=allnodes\n');
%     fprintf(u3,'U\n');
%     
% fprintf(u3,'*END STEP\n');

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
    
    if k==1
         fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
    else
         fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=MOD\n');
    end
    fprintf(u2,'Part-1-1.%d, PY, %f\n',[Pd(:,1), lambda(k)*Pd(:,2)]');
        
    fprintf(u3,'** \n');
    fprintf(u3,'** OUTPUT REQUESTS\n');
    fprintf(u3,'** \n');
    fprintf(u3,'*Restart, write, frequency=0\n');
    fprintf(u3,'** \n');
    fprintf(u3,'** FIELD OUTPUT: F-Output-1\n');
    fprintf(u3,'** \n');
    fprintf(u3,'*Output, field, variable=PRESELECT\n');
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
    
    
% fprintf(u3,'** ----------------------------------------------------------------\n');
% fprintf(u3,['*Step, name=Eig-',num2str(stepnum),', nlgeom=NO, perturbation\n']);
% fprintf(u3,'*Buckle, eigensolver=lanczos\n');
% fprintf(u3,'10, , , , \n');
% 
%     
%     fprintf(u3,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
% %     fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), (1-lambda(k))*Pd(:,2)]');
%     fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), Pd(:,2)]');
%     fprintf(u3,'*NODE PRINT, nset=allnodes\n');
%     fprintf(u3,'U\n');
% fprintf(u3,'*END STEP\n');
    
end

fclose(u3);

% u4 = fopen(['AbaqusRuns/',filename,'_initKg.inp'],'w');
% 
% fprintf(u4,['*Include, input=',filename,'-model.inp\n']);
% fprintf(u4,'** ----------------------------------------------------------------\n');
% fprintf(u4,'*STEP, name=Lambda-1\n');
% fprintf(u4,'*MATRIX GENERATE, STIFFNESS\n');
% fprintf(u4,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
% fprintf(u4,'*END STEP\n');
% 
%     stepnum = 1;
%     fprintf(u4,'** ----------------------------------------------------------------\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,['** STEP: Step-',num2str(stepnum),'\n']);
%     fprintf(u4,'** \n');
%     fprintf(u4,['*Step, name=Step-',num2str(stepnum),', nlgeom=YES\n']);
%     fprintf(u4,'*Static\n');
%     fprintf(u4,['1',', ','1',', ','0.0001',', ','1','\n']);
% 
%     fprintf(u4,'** \n');
%     fprintf(u4,'** LOADS\n');
%     fprintf(u4,'** \n');
% 
%     fprintf(u4,'*Dload, CONSTANT RESULTANT=YES, OP=NEW\n');
%     fprintf(u4,'Part-1-1.%d, PY, %f\n',[Pd(:,1), Pd(:,2)/abs(p/1000)]');
%         
%     fprintf(u4,'** \n');
%     fprintf(u4,'** OUTPUT REQUESTS\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,'*Restart, write, frequency=0\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,'** FIELD OUTPUT: F-Output-1\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,'*Output, field, variable=PRESELECT\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,'** HISTORY OUTPUT: H-Output-1\n');
%     fprintf(u4,'** \n');
%     fprintf(u4,'*Output, history, variable=PRESELECT\n');
%     fprintf(u4,'*NODE PRINT, nset=allnodes\n');
%     fprintf(u4,'U\n');
%     fprintf(u4,'*EL PRINT\n');
%     fprintf(u4,'\n');
%     fprintf(u4,'SE\n');
%     fprintf(u4,'*EL PRINT\n');
%     fprintf(u4,'\n');
%     fprintf(u4,'SF\n');
%     fprintf(u4,'*End Step\n');
%     
%     stepnum = stepnum + 1;
%     fprintf(u4,'** ----------------------------------------------------------------\n');
%     fprintf(u4,['*STEP, name=Lambda-',num2str(stepnum),'\n']);
%     fprintf(u4,'*MATRIX GENERATE, STIFFNESS\n');
%     fprintf(u4,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
%     fprintf(u4,'*END STEP\n');
% 
% fclose(u4);

end


function plotMesh(Nodes,Elements)
   hold on
   for i = 1:size(Elements,1)
       nnums = Elements(i,2:end);
       nnums = [nnums,nnums(1)];
       xcoords = Nodes(nnums,2);
       ycoords = Nodes(nnums,3);
       zcoords = Nodes(nnums,4);
       plot3(xcoords,ycoords,zcoords,'bx-')
   end
end