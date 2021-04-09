function [filename,lambda,BC,Nodes,Elements,P,dofpNode]  = detKt2D(L,numofelm,lambda,loadFactor,eltype,~,modelprops,AbaqusRunsFolder)
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
  eltype = 'B32OSH';
 end

 
 
 if lambda(1) == 0
  lambda(1) = [];
 end
 
 % pure SI units: Newtons, meters, Pascals, etc.
 filename = ['detKt-',eltype,'-',num2str(numofelm(1)),'-',num2str(numofelm(end)),'-len-',num2str(L),'-f-',num2str(loadFactor),'-eps',num2str(modelprops.epsilon)];
 
 %% IPE400
 h = (400)*10^(-3); %[m]
 b = 180e-3; %[m]
 tw = 8.6e-3; %[m]
 tf = 13.5e-3; %[m]
 

 
 %Iy = 2*(tf^3*b/12 + tf*b*(h/2)^2) + (h)^3*tw/12;
 
 %% Load
 P = loadFactor*1e6; %[N?]
 %M = 0; % [N m?]
 
%% Finite Element Model
 
 xcoords = unique([linspace(-L,0,numofelm(1)+1)';linspace(0,L,numofelm(end)+1)']);
 ycoords = 0*xcoords;
 %zcoords = 0*xcoords;
 xcoords(abs(xcoords)<1e-12) = 0;
 
 Nodes = [ctranspose(1:length(xcoords)), xcoords, ycoords];
 Elements = [ctranspose(1:size(Nodes,1)-1),Nodes(1:end-1,1),Nodes(2:end,1)];
 
 rpLeft = Nodes(1,1);
 rpMiddle= Nodes(Nodes(:,2)==0,1);
 rpRight = Nodes(end,1);
 
 if strcmpi('B32',eltype(1:3)) || strcmpi('B22',eltype(1:3))
  Nodes2 = 0.5*(Nodes(1:end-1,2:end) + Nodes(2:end,2:end));
  Nodes2 = [Nodes(end,1) + ctranspose(1:size(Nodes2,1)),Nodes2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];%Nr,Anfang,mitte,Ende
 end

   midnode1 = Nodes(Nodes(:,2)==0,1);
%   Nodes = [Nodes; Nodes(midnode1,:)];
%   %midnode2 = size(Nodes,1);
   midnode2=rpMiddle;
%   Nodes(end,1) = midnode2;
   Elements(Elements(:,2)==midnode1,2) = midnode2;
 
 %% Boundary conditions
 if strcmp(eltype,'B21') || strcmp(eltype,'B23') || strcmp(eltype,'B22')
  dofpNode=6;
 else
  error('MyProgram:Element','unknown Element')
 end
 %+1...Nx
 %+2...Ny
 %+6...Mz
 BCrpLeft=[1,2];
 BCrpMiddle=[2,6];
 BCmidnode2=BCrpMiddle;
 BCrpRight=[1,2];
 BC = ...
 [dofpNode*(rpLeft  - 1) + 1, 0;
  dofpNode*(rpLeft  - 1) + 2, 0;
  %dofpNode*(rpLeft  - 1) + 3, 0;
  %dofpNode*(rpLeft  - 1) + 4, 0;
  %dofpNode*(rpMiddle- 1) + 1, 0;
  dofpNode*(rpMiddle- 1) + 2, 0;
  %dofpNode*(rpMiddle- 1) + 3, 0;
  %dofpNode*(rpMiddle- 1) + 4, 0;
  dofpNode*(rpMiddle- 1) + 6, 0;
  %dofpNode*(midnode2- 1) + 1, 0;
  %dofpNode*(midnode2- 1) + 2, 0;
  %dofpNode*(midnode2- 1) + 3, 0;
  %dofpNode*(midnode2- 1) + 4, 0;
  %dofpNode*(rpRight - 1) + 1, 0;
  dofpNode*(rpRight - 1) + 2, 0;
  %dofpNode*(rpRight - 1) + 3, 0;
  %dofpNode*(rpRight - 1) + 4, 0
  ];
 for i=1:numel(BCrpLeft)
  BC=[BC;dofpNode*(rpLeft - 1) + BCrpLeft(i), 0]; %#ok<AGROW>
 end
 for i=1:numel(BCrpMiddle)
  BC=[BC;dofpNode*(rpMiddle - 1) + BCrpMiddle(i), 0]; %#ok<AGROW>
 end
 for i=1:numel(BCmidnode2)
  BC=[BC;dofpNode*(midnode2 - 1) + BCmidnode2(i), 0]; %#ok<AGROW>
 end
 for i=1:numel(BCrpRight)
  BC=[BC;dofpNode*(rpRight - 1) + BCrpRight(i), 0]; %#ok<AGROW>
 end
 [~,idx]=unique(BC(:,1));
 BC=BC(idx,:);
 
 u1 = fopen([AbaqusRunsFolder,filename,'.inpData'],'w');
 if u1==-1
  error('MyProgram:FileNotOpen','kann die Datei nicht oeffnen')
 else
  fprintf(u1,'*Part, name=Part-1\n');
  fprintf(u1,'*Node\n');
  if strcmpi(eltype(1:2),'B2')
   fprintf(u1,'%d, %f, %f\n',Nodes');
  else
   fprintf(u1,'%d, %f, %f, %f\n',Nodes');
  end
  fprintf(u1,['*Element, type=',eltype,'\n']);
  if strcmpi(eltype(1:3),'B32') || strcmpi(eltype(1:3),'B22')
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
  fprintf(u1,'*Nset, nset=Mitte, instance=Part-1-1\n');
  fprintf(u1,[num2str(rpMiddle),', ',num2str(midnode2),'\n']);
  fprintf(u1,'*Nset, nset=MomemntMitte, instance=Part-1-1\n');
  fprintf(u1,[num2str(midnode2),'\n']);
  fprintf(u1,'*Nset, nset=ForceMitte, instance=Part-1-1\n');
  fprintf(u1,[num2str(midnode2),'\n']);
  fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
  fprintf(u1,[num2str(rpRight),'\n']);
  
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+1,midnode1,midnode2]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');
  
  fprintf(u1,'*End Assembly\n');
  
  fprintf(u1,'*Material, name=Material-1\n');
  fprintf(u1,'*Elastic\n');
  fprintf(u1,'2.1e+11, 0.3\n');
  
  %% Boundary conditions
  %%
  
  fprintf(u1,'*Boundary\n');
%   %fprintf(u1,'leftend,  1, 1\n');
%   %fprintf(u1,'leftend,  2, 2\n');
%   fprintf(u1,'leftend,  3, 3\n');
%   fprintf(u1,'leftend,  4, 4\n');
%   %fprintf(u1,'Mitte,  1, 1\n');
%   %fprintf(u1,'Mitte,  2, 2\n');
%   fprintf(u1,'Mitte,  3, 3\n');
%   fprintf(u1,'Mitte,  4, 4\n');
%   %fprintf(u1,'Mitte,  6, 6\n');
%   %fprintf(u1,'rightend,  1, 1\n');
%   %fprintf(u1,'rightend,  2, 2\n');
%   fprintf(u1,'rightend,  3, 3\n');
%   fprintf(u1,'rightend,  4, 4\n');
 for i=1:numel(BCrpLeft)
  FoN=BCrpLeft(i);
  fprintf(u1,'leftend,  %d, %d\n',FoN,FoN);
 end
 for i=1:numel(BCrpMiddle)
  FoN=BCrpMiddle(i);
  fprintf(u1,'Mitte,  %d, %d\n',FoN,FoN);
 end
 for i=1:numel(BCrpRight)
  FoN=BCrpRight(i);
  fprintf(u1,'rightend,  %d, %d\n',FoN,FoN);
 end
  
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
   %fprintf(u3,'leftend, 1, %f\n',lambda(k)*P);
   fprintf(u3,'ForceMitte, 1, %f\n',-lambda(k)*P);
   %fprintf(u3,'MomemntMitte, 5, %f\n',lambda(k)*M);
   %fprintf(u3,'rightend, 5, %f\n',-lambda(k)*M);
   
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
