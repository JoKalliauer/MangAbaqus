function [filename,lambda,BC,Nodes,Elements,load,dofpNode]  = eccenCompressionBeam(L,numofelm,lambda,loadFactor,elType,ecc,modelprops,AbaqusRunsFolder)
 if nargin<1
  L = 5.0;
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
 if lambda(1) == 0
  lambda(1) = [];
 end
 
 % pure SI units: Newtons, meters, Pascals, etc.
 filename =  ['ecc-',elType,'-',num2str(numofelm(end)),'-l',num2str(L),'-e',num2str(ecc),'-f',num2str(loadFactor),'-eps',num2str(modelprops.epsilon)];
 
 %% IPE400
 h = (400)*10^(-3); %[m]
 b = 180e-3; %[m]
 tw = 8.6e-3; %[m]
 tf = 13.5e-3; %[m]
 
 if nargin<6
  ecc = 0.05; % eccectricity
 else
  assert(numel(ecc)==1,'dimension of ecc must be one');
 end
 
 %Iy = 2*(tf^3*b/12 + tf*b*(h/2)^2) + (h)^3*tw/12;
 
 %% Load
 P = loadFactor*500e3; %[N?]
% <<<<<<< HEAD
 load=lambda*P;
% =======
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
 %% Finite Element Model
 
 xcoords = linspace(-L/2,L/2,numofelm+1)';
 ycoords = 0*xcoords;
 zcoords = 0*xcoords;
 xcoords(abs(xcoords)<1e-12) = 0;
 
% <<<<<<< HEAD
 if ecc~=0
% =======
%  %if ecc~=0
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 % add eccentric rods:
 xcoords = [xcoords(1); xcoords; xcoords(end)];
 ycoords = [ecc; ycoords; ecc];
 zcoords = [0; zcoords; 0];
% <<<<<<< HEAD
 end
 
 Nodes = [ctranspose(1:length(xcoords)), xcoords, ycoords, zcoords];
 Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
 

 
 if ecc~=0
  rpLeft1 = Nodes(2,1);
  rpRight1 = Nodes(end-1,1);
  rpLeft2 = Nodes(1,1); %LEFTENDFORCE
  rpRight2 = Nodes(end,1);
 else
  rpLeft1 = Nodes(1,1);
  rpRight1 = Nodes(end,1);
% =======
%  %end
%  
%  Nodes = [ctranspose(1:length(xcoords)), xcoords, ycoords, zcoords];
%  Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
%  
%  rpLeft1 = Nodes(2,1);
%  rpRight1 = Nodes(end-1,1);
%  
%  if ecc~=0
%   rpLeft2 = Nodes(1,1); %LEFTENDFORCE
%   rpRight2 = Nodes(end,1);
%  else
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  rpLeft2 = rpLeft1; %LEFTENDFORCE=leftend
  rpRight2 = rpRight1;
 end
 
 if strcmpi('B32',elType(1:3))
  Nodes2 = 0.5*(Nodes(1:end-1,2:end) + Nodes(2:end,2:end));
  Nodes2 = [Nodes(end,1) + ctranspose(1:size(Nodes2,1)),Nodes2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 end
 
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
 fprintf(u1,'%d, %f, %f, %f\n',Nodes');
 fprintf(u1,['*Element, type=',elType,'\n']);
% <<<<<<< HEAD
 if ecc~=0
  if strcmpi(elType(1:3),'B32')
   fprintf(u1,'%d, %d, %d, %d\n',Elements(2:end-1,:)');
  else
   fprintf(u1,'%d, %d, %d\n',Elements(2:end-1,:)');
  end
% =======
%  if strcmpi(elType(1:3),'B32')
%   fprintf(u1,'%d, %d, %d, %d\n',Elements(2:end-1,:)');
%  else
%   fprintf(u1,'%d, %d, %d\n',Elements(2:end-1,:)');
%  end
%  if ecc~=0
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
  fprintf(u1,['*Element, type=',elType,'\n']);
  if strcmpi(elType(1:3),'B32')
   fprintf(u1,'%d, %d, %d, %d\n',[Elements(1,:); Elements(end,:)]');
  else
   fprintf(u1,'%d, %d, %d\n',[Elements(1,:); Elements(end,:)]');
  end
% <<<<<<< HEAD
 else
  if strcmpi(elType(1:3),'B32')
   fprintf(u1,'%d, %d, %d, %d\n',Elements(1:end,:)');
  else
   fprintf(u1,'%d, %d, %d\n',Elements(1:end,:)');
  end
 end
 
 fprintf(u1,'*Elset, elset=AllElements\n');
 if ecc~=0
  fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(2:end-1,1));
 else
  fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(1:end,1));
 end
 if length(Elements(:,1))/16~=floor(length(Elements(:,1))/16)
  fprintf(u1,'\n');
 end
 if ecc~=0
  fprintf(u1,'*Elset, elset=Vert\n');
  fprintf(u1,'%d, %d\n',[Elements(1,1),Elements(end,1)]);
 end
% =======
%  end
%  
%  fprintf(u1,'*Elset, elset=AllElements\n');
%  fprintf(u1,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',Elements(2:end-1,1));
%  if length(Elements(:,1))/16~=floor(length(Elements(:,1))/16)
%   fprintf(u1,'\n');
%  end
%  fprintf(u1,'*Elset, elset=Vert\n');
%  fprintf(u1,'%d, %d\n',[Elements(1,1),Elements(end,1)]);
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 
 fprintf(u1,'** Section: Section-1  Profile: Profile-1\n');
 fprintf(u1,'*Beam Section, elset=AllElements, material=Material-1, temperature=VALUES, section=I\n');
 profile = [h/2, h, b, b, tf, tf, tw];
 fprintf(u1,'%f, %f, %f, %f, %f, %f, %f \n',profile);
 fprintf(u1,'0.,-1.,0.\n');
 fprintf(u1,'** Section: Section-2  Profile: Profile-2\n');
 fprintf(u1,'*Beam Section, elset=Vert, material=Material-1, temperature=VALUES, section=I\n');
 fprintf(u1,'1.0, 2.0, 1.0, 1.0, 0.1, 0.1, 0.1 \n');
 fprintf(u1,'1.,0.,0.\n');
 fprintf(u1,'*End Part\n');
 
 fprintf(u1,'*Assembly, name=Assembly\n');
 fprintf(u1,'*Instance, name=Part-1-1, part=Part-1\n');
 fprintf(u1,'*End Instance\n');
 
 fprintf(u1,'*Nset, nset=allnodes, instance=Part-1-1, generate\n');
 fprintf(u1,['1,',num2str(Nodes(end,1)),',1\n']);
 
 fprintf(u1,'*Nset, nset=leftend, instance=Part-1-1\n');
 fprintf(u1,[num2str(rpLeft1),'\n']);
 fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
 fprintf(u1,[num2str(rpRight1),'\n']);
 
 fprintf(u1,'*Nset, nset=LEFTENDFORCE, instance=Part-1-1\n');
 fprintf(u1,[num2str(rpLeft2),'\n']);
 fprintf(u1,'*Nset, nset=rightendforce, instance=Part-1-1\n');
 fprintf(u1,[num2str(rpRight2),'\n']);
 
 fprintf(u1,'*End Assembly\n');
 
 fprintf(u1,'*Material, name=Material-1\n');
 fprintf(u1,'*Elastic\n');
% <<<<<<< HEAD
 fprintf(u1,'2.1e+11, 0.3\n'); %[N/m² , -]?
 
 %% Boundary conditions
 %dofpNode=7;
 if strcmp(elType,'B32') || strcmp(elType,'B32H') || strcmp(elType,'B31') || strcmp(elType,'B33') ||  strcmp(elType,'B31H') || strcmp(elType,'B33H')
  dofpNode=6;
% =======
%  fprintf(u1,'2.1e+11, 0.3\n');
%  
%  %% Boundary conditions
%  dofpNode=7;
%  if strcmp(elType,'B32OS')
%   %dofpNode=6;
% >>>>>>> c8d007979d050d2fdcd2c9ed43fa8f6b3bcff9d2
 else
  dofpNode=7;
 end
 BC = [dofpNode*(rpLeft1 - 1) + 1, 0
       dofpNode*(rpLeft1 - 1) + 2, 0;
       dofpNode*(rpLeft1 - 1) + 3, 0;
       dofpNode*(rpLeft1 - 1) + 4, 0;
       dofpNode*(rpRight1 - 1) + 2, 0;
       dofpNode*(rpRight1 - 1) + 3, 0;
       dofpNode*(rpRight1 - 1) + 4, 0];
 %%
 
 fprintf(u1,'*Boundary\n');
 fprintf(u1,'leftend,  1, 1\n');
 fprintf(u1,'leftend,  2, 2\n');
 fprintf(u1,'leftend,  3, 3\n');
 fprintf(u1,'leftend,  4, 4\n');
 fprintf(u1,'rightend,  2, 2\n');
 fprintf(u1,'rightend,  3, 3\n');
 fprintf(u1,'rightend,  4, 4\n');
 
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
  fprintf(u3,'LEFTENDFORCE, 1, %f\n',lambda(k)*P);
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
