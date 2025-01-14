function [filename,lambda,BC,Nodes,Elements]  = TL_arch(~,numofelm,lambda,RefLast,elType,AbaqusRunsFolder,modelprops)
% if nargin<1
%     dummy = [];
% end
if nargin<2
    numofelm = 20;
end
if nargin<3
    lambda = 0:0.1:1;
end

if nargin<5
    elType = 'B22H';
end
if lambda(1) == 0
    lambda(1) = [];
end

%% RECT
 h = 20e-2; %[m]
 b = 10e-2; %[m]
  
%% Span
 L = 600e-2; % m
 H = 240e-2; % m


% pure SI units: Newtons, meters, Pascals, etc.
filename = ['TL_arch-',elType,'-',num2str(numofelm(end)),'-RefLast-',num2str(RefLast),'-eps',num2str(modelprops.epsilon)];

%% Load
 p = -RefLast; %[N/m]
 
 %% Finite Elements Size
 resolution = L/double(numofelm(end));
 
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
    
 Nodes = [ctranspose(1:length(xcoords)), xcoords, ycoords];
 Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
 
 rpLeft = Nodes(1,1);
%  leftnodes = Nodes(1,1);
 rpRight = Nodes(end,1);
%  rightnodes = Nodes(end,1);
 
 if strcmpi('B22',elType(1:3))  
    Nodes2 = 0.5*(Nodes(1:end-1,2:end) + Nodes(2:end,2:end));
    Nodes2 = [Nodes(end,1) + ctranspose(1:size(Nodes2,1)),Nodes2];
    Nodes = [Nodes; Nodes2];
    Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 end
 
 P = zeros(size(Nodes,1),1);
 Pd = zeros(size(Elements,1),1);
 for i = 1:size(Elements,1)
     node1 = Elements(i,2); node2 = Elements(i,3);
     if strcmpi(elType(1:3),'B22')
         node3 = Elements(i,4);
         coords = [Nodes(node1,2:3); Nodes(node2,2:3); Nodes(node3,2:3)];
            % extract position on XY plane
           x1 = coords(1,1); x2 = coords(2,1); x3 = coords(3,1); 
           y1 = coords(1,2);
           %y2 = coords(2,2);
           y3 = coords(3,2);

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
 
 %P = [Nodes(:,1),P];
 Pd = [Elements(:,1),Pd];
 
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
fprintf(u1,'%d, %f, %f \n',Nodes');
fprintf(u1,['*Element, type=',elType,'\n']);
if strcmpi(elType(1:3),'B22')
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

fprintf(u1,'*End Assembly\n');

fprintf(u1,'*Material, name=Material-1\n');
fprintf(u1,'*Elastic\n');
fprintf(u1,'2e+11, 0.3\n');

%% Boundary conditions
 BC = [6*(rpLeft - 1) + 1, 0
       6*(rpLeft - 1) + 2, 0;
       6*(rpRight - 1) + 1, 0;
       6*(rpRight - 1) + 2, 0];
%%

fprintf(u1,'*Boundary\n');
fprintf(u1,'leftend,  1, 1\n');
fprintf(u1,'leftend,  2, 2\n');
fprintf(u1,'rightend,  1, 1\n');
fprintf(u1,'rightend,  2, 2\n');

u3 = fopen([AbaqusRunsFolder,filename,'.inp'],'w');

fprintf(u3,['*Include, input=',filename,'.inpData\n']);
fprintf(u3,'** ----------------------------------------------------------------\n');
fprintf(u3,'*STEP, name=Lambda-1\n');
fprintf(u3,'*MATRIX GENERATE, STIFFNESS\n');
fprintf(u3,'*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n');
fprintf(u3,'*END STEP\n');

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
    fprintf(u3,'Part-1-1.%d, PY, %f\n',[Pd(:,1), lambda(k)*Pd(:,2)]');
        
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
   
end

fclose(u3);

end