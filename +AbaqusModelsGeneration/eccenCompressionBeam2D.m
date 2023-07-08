function [filename,lambda,BC,Nodes,Elements,load,dofpNode,h]  = eccenCompressionBeam2D(L0,numofelm,lambda,loadFactor,elType,ecc,modelprops,AbaqusRunsFolder)
if nargin<1
 L0 = 5.0;
end
if nargin<2
 numofelm = 20;
end
if nargin<3
 lambda = 0:0.1:1;
end
if nargin<4
 loadFactor = 1.0;
end
if nargin<5
 elType = 'B32OSH';
end

MV=modelprops.MeterValue;
assert(MV==1,'MV different to one not tested')

if lambda(1) == 0
 lambda(1) = [];
end

% pure SI units: Newtons, meters, Pascals, etc.
filename = ['eccenCompressionBeam2D-',elType,'-',num2str(numofelm(end)),'-len-',num2str(L0),'-ecc-',num2str(ecc),'-loadfac-',num2str(loadFactor)];

L=L0*MV;
%% IPE400
h = (400)*10^(-3)*MV;
b = 180e-3*MV;
tw = 8.6e-3*MV;
tf = 13.5e-3*MV;

if nargin<6
 ecc = 0.05; % eccectricity
end
eccMV=ecc*MV;

%Iy = 2*(tf^3*b/12 + tf*b*(h/2)^2) + (h)^3*tw/12;

%% Load
P = loadFactor*1.360246541947329e+07;
load=lambda*P;

%% Finite Element Model

xcoords = linspace(-L/2,L/2,numofelm+1)';
ycoords = 0*xcoords;
zcoords = 0*xcoords;
xcoords(abs(xcoords)<1e-12) = 0;

if eccMV~=0
 % add eccentric rods:
 xcoords = [xcoords(1); xcoords; xcoords(end)];
 ycoords = [eccMV; ycoords; eccMV];
 zcoords = [0; zcoords; 0];
end

Nodes1 = [ctranspose(1:length(xcoords)), xcoords, ycoords, zcoords];
Elements = [ctranspose(1:size(Nodes1,1)-1),Nodes1(1:end-1,1),Nodes1(2:end,1)];

if eccMV~=0
 rpLeft1 = Nodes1(2,1);
 rpRight1 = Nodes1(end-1,1);
 rpLeft2 = Nodes1(1,1);
 rpRight2 = Nodes1(end,1);
else
 rpLeft1 = Nodes1(1,1);
 rpRight1 = Nodes1(end,1);
 rpLeft2 = rpLeft1; %LEFTENDFORCE=leftend
 rpRight2 = rpRight1;
end

if strcmpi('B22',elType(1:3))
 Nodes2 = 0.5*(Nodes1(1:end-1,2:end) + Nodes1(2:end,2:end));
 Nodes2 = [Nodes1(end,1) + ctranspose(1:size(Nodes2coord,1)),Nodes2];
 Nodes = [Nodes1; Nodes2];
 Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 if eccMV~=0
  OrderedBeamNodes(1:2:2*numofelm+1)=Nodes1(2:end-1,1);
  OrderedBeamNodes(2:2:2*numofelm)=Nodes2(2:end-1,1);
 else
  OrderedBeamNodes=Nodes1(1:end,1);
 end
else
 Nodes=Nodes1;
 if eccMV~=0
  OrderedBeamNodes=Nodes1(2:end-1,1);
 else
  OrderedBeamNodes=Nodes1(1:end,1);
 end
end

u1 = fopen([AbaqusRunsFolder,filename,'.inpData'],'w');

fprintf(u1,'*Part, name=Part-1\n');
fprintf(u1,'*Node\n');
fprintf(u1,'%d, %f, %f, %f\n',Nodes');
fprintf(u1,['*Element, type=',elType,'\n']);%Stab-Elemente
if eccMV~=0
 if strcmpi(elType(1:3),'B32')
  fprintf(u1,'%d, %d, %d, %d\n',Elements(2:end-1,:)');
 else
  fprintf(u1,'%d, %d, %d\n',Elements(2:end-1,:)');
 end
 fprintf(u1,['*Element, type=',elType,'\n']);%belastungselemente
 if strcmpi(elType(1:3),'B32')
  fprintf(u1,'%d, %d, %d, %d\n',[Elements(1,:); Elements(end,:)]');
 else
  fprintf(u1,'%d, %d, %d\n',[Elements(1,:); Elements(end,:)]');
 end
else
 if strcmpi(elType(1:3),'B32')
  fprintf(u1,'%d, %d, %d, %d\n',Elements(1:end,:)');
 else
  fprintf(u1,'%d, %d, %d\n',Elements(1:end,:)');
 end
end

fprintf(u1,'*Elset, elset=AllElements\n');
if eccMV~=0
 fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(2:end-1,1));
else
 fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(1:end,1));
end
if length(Elements(:,1))/16~=floor(length(Elements(:,1))/16)
 fprintf(u1,'\n');
end
if eccMV~=0
 fprintf(u1,'*Elset, elset=Vert\n');
 fprintf(u1,'%d, %d\n',[Elements(1,1),Elements(end,1)]);
end

fprintf(u1,'** Section: Section-1  Profile: Profile-1\n');
fprintf(u1,'*Beam Section, elset=AllElements, material=Material-1, temperature=VALUES, section=I\n');
profile = [h/2, h, b, b, tf, tf, tw];
fprintf(u1,'%f, %f, %f, %f, %f, %f, %f \n',profile);
fprintf(u1,'0.,0.,-1.\n');
fprintf(u1,'** Section: Section-2  Profile: Profile-2\n');
fprintf(u1,'*Beam Section, elset=Vert, material=Material-1, temperature=VALUES, section=I\n');
fprintf(u1,'1.0, 2.0, 1.0, 1.0, 0.1, 0.1, 0.1 \n');
fprintf(u1,'1.,0.,0.\n');
fprintf(u1,'*End Part\n');

fprintf(u1,'*Assembly, name=Assembly\n');
fprintf(u1,'*Instance, name=Part-1-1, part=Part-1\n');
fprintf(u1,'*End Instance\n');

fprintf(u1,'*Nset, nset=allnodes, instance=Part-1-1, generate\n');
fprintf(u1,['1,',num2str(max(Nodes(:,1))),',1\n']);
fprintf(u1,'*Nset, nset=OrderedBeamNodes, instance=Part-1-1\n');
fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',OrderedBeamNodes);%per line only 16 Values read by abaqus
%for n=20: 2,25,3,26,4,27,5,28,6,29,7,30,8,31,9,32,10,33,11,34,12,35,13,36,14,37,15,38,16,39,17,40,18,41,19,42,20,43,21,44,22
fprintf(u1,'\n');%can lead to two linebreaks but I dont care

fprintf(u1,'*Nset, nset=leftend, instance=Part-1-1\n');
fprintf(u1,[num2str(rpLeft1),'\n']);
fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
fprintf(u1,[num2str(rpRight1),'\n']);

fprintf(u1,'*Nset, nset=leftendforce, instance=Part-1-1\n');
fprintf(u1,[num2str(rpLeft2),'\n']);
fprintf(u1,'*Nset, nset=rightendforce, instance=Part-1-1\n');
fprintf(u1,[num2str(rpRight2),'\n']);

fprintf(u1,'*End Assembly\n');

fprintf(u1,'*Material, name=Material-1\n');
fprintf(u1,'*Elastic\n');
fprintf(u1,'2.1e+11, 0.3\n');

%% Boundary conditions
if strcmp(elType,'B21H')
 dofpNode=6;
elseif strcmp(elType,'B22OS') || strcmp(elType,'B22OSH') || strcmp(elType,'B21OS') || strcmp(elType,'B21OSH')
 dofpNode=7;
elseif strcmp(elType,'B22') || strcmp(elType,'B22H') || strcmp(elType,'B21') || strcmp(elType,'B23') || strcmp(elType,'B23H')
 %dofpNode=6;
 error('MyPrgm:Element','element not tested')
else
 %dofpNode=6;
 error('MyPrgm:Element','unknown element')
end
BC = [dofpNode*(rpLeft1 - 1) + 1, 0
 dofpNode*(rpLeft1 - 1) + 2, 0;
 dofpNode*(rpRight1 - 1) + 2, 0;
 %dofpNode*(Nodes(:,1)-1) + 3, zeros(size(Nodes,1),1);
 %dofpNode*(Nodes(:,1)-1) + 4, zeros(size(Nodes,1),1);
 %dofpNode*(Nodes(:,1)-1) + 5, zeros(size(Nodes,1),1);
 %dofpNode*(Nodes(:,1)-1) + 7, zeros(size(Nodes,1),1)
 ];
%%

fprintf(u1,'*Boundary\n');
fprintf(u1,'leftend,  1, 1\n');
fprintf(u1,'leftend,  2, 2\n');
fprintf(u1,'rightend,  2, 2\n');
fprintf(u1,'allnodes,  3, 3\n');
fprintf(u1,'allnodes,  4, 4\n');
fprintf(u1,'allnodes,  5, 5\n');
fprintf(u1,'allnodes,  7, 7\n');

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
  fprintf(u3,'*Cload, OP=NEW\n');
 else
  fprintf(u3,'*Cload, OP=MOD\n');
 end
 fprintf(u3,'leftendforce, 1, %f\n',lambda(k)*P);
 fprintf(u3,'rightendforce, 1, %f\n',-lambda(k)*P);
 
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
