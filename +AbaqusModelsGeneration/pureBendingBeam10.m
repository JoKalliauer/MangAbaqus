function [filename,lambda,BC,Nodes,Elements,Last,dofpNode,h]  = pureBendingBeam10(Linp,numofelm,lambda,loadFactor,elType,modelprops,AbaqusRunsFolder)
 if nargin<1
  Linp = 5.0;
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
 if sum(strcmp(fieldnames(modelprops), 'orientate')) == 0
  modelprops.orientate=5;
 end
 assert(ismember(modelprops.orientate,[5 6]),'modelprops.orientate must be 5=y or 6=z')
 
 % pure SI units: Newtons, meters, Pascals, etc.
 
  MV=modelprops.MeterValue;
  

  

 filename = ['B10B',num2str(modelprops.orientate),'-',elType,num2str(numofelm(1)),'-l',num2str(Linp),'-f',num2str(loadFactor),'-eps',num2str(modelprops.epsilon),'-u',num2str(MV)];

 
 %% IPE400
 h = (400)*10^(-3)*MV; %[m]
 b = 180e-3*MV; %[m]
 tw = 8.6e-3*MV; %[m]
 tf = 13.5e-3*MV; %[m]
 L=Linp*MV; %Laenge in [m]
 
 %Iy = 2*(tf^3*b/12 + tf*b*(h/2)^2) + (h)^3*tw/12;
 
 Emodul=2.1e+11/MV;
 
 %% Load
 M = loadFactor*5e6*MV*MV; %[N*m ?]
 Last=lambda*M;
 
 %% Finite Element Model
 
 xcoords = linspace(-L/2,L/2,numofelm(1)+1)';
 ycoords = 0*xcoords;
 zcoords = 0*xcoords;
 xcoords(abs(xcoords)<1e-12) = 0;
 
 Nodes = [ctranspose(1:length(xcoords)), xcoords, ycoords, zcoords];
 Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
 
 rpLeft = Nodes(1,1);
 rpRight = Nodes(end,1);
 
 if strcmpi('B32',elType(1:3))
  Nodes2 = 0.5*(Nodes(1:end-1,2:end) + Nodes(2:end,2:end));
  Nodes2 = [Nodes(end,1) + ctranspose(1:size(Nodes2,1)),Nodes2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];
 end
 
 u1 = fopen([AbaqusRunsFolder,filename,'-model.inp'],'w');
 if u1==-1
  warning('MyProgram:FileNotOpen','kann die Datei nicht oeffnen')
 else
  fprintf(u1,'*Part, name=Part-1\n');
  fprintf(u1,'*Node\n');
  fprintf(u1,'%d, %f, %f, %f\n',Nodes');
  fprintf(u1,['*Element, type=',elType,'\n']);
  if strcmpi(elType(1:3),'B32')
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
  fprintf(u1,'*Beam Section, elset=AllElements, material=Material-1, temperature=VALUES, section=I\n');
  profile = [h/2, h, b, b, tf, tf, tw];
  fprintf(u1,'%f, %f, %f, %f, %f, %f, %f \n',profile);
  fprintf(u1,'0.,-1.,0.\n');
  fprintf(u1,'*End Part\n');
  
  fprintf(u1,'*Assembly, name=Assembly\n');
  fprintf(u1,'*Instance, name=Part-1-1, part=Part-1\n');
  fprintf(u1,'*End Instance\n');
  
  fprintf(u1,'*Nset, nset=allnodes, instance=Part-1-1, generate\n');
  fprintf(u1,['1,',num2str(Nodes(end,1)),',1\n']);
  
  fprintf(u1,'*Nset, nset=leftend, instance=Part-1-1\n');
  fprintf(u1,[num2str(rpLeft),'\n']);
  fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
  fprintf(u1,[num2str(rpRight),'\n']);
  
  fprintf(u1,'*End Assembly\n');
  
  fprintf(u1,'*Material, name=Material-1\n');
  fprintf(u1,'*Elastic\n');
%   fprintf(u1,'2.1e+11, 0.3\n');
  fprintf(u1,[num2str(Emodul) ', 0.3\n']);
  
  %% Boundary conditions
  if strcmp(elType,'B32') || strcmp(elType,'B31') || strcmp(elType,'B33') ||  strcmp(elType,'B31H') || strcmp(elType,'B33H') || strcmp(elType,'B32H')
   dofpNode=6;
  elseif strcmp(elType,'B32OS') || strcmp(elType,'B32OSH')|| strcmp(elType,'B31OS') || strcmp(elType,'B31OSH')
   dofpNode=7;
  else
   %dofpNode=7;
   dofpNode=[];
  end
  BC = [dofpNode*(rpLeft - 1) + 1, 0
        dofpNode*(rpLeft - 1) + 2, 0;
        dofpNode*(rpLeft - 1) + 3, 0;
        dofpNode*(rpLeft - 1) + 4, 0;
        dofpNode*(rpRight - 1) + 2, 0;
        dofpNode*(rpRight - 1) + 3, 0;
        dofpNode*(rpRight - 1) + 4, 0];
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
  
  fprintf(u3,['*Include, input=',filename,'-model.inp\n']);
  
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
   fprintf(u3,'leftend, %d, %f\n',modelprops.orientate,lambda(k)*M);
   fprintf(u3,'rightend, %d, %f\n',modelprops.orientate,-lambda(k)*M);
   
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
end
