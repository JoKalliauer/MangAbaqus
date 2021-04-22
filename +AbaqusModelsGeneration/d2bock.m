function [filename,lambda,BC,Nodes,Elements,P,dofpNode]  = d2bock(L0,numofelm,lambda,loadFactor,eltype,modelprops,AbaqusRunsFolder)
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
  eltype = 'B32OSH';
 end

 MV=modelprops.MeterValue;
 
 if lambda(1) == 0
  lambda(1) = [];
 end
 
 % pure SI units: Newtons, meters, Pascals, etc.
 filename = ['d2bock-',eltype,'-',num2str(numofelm(1)),'-',num2str(numofelm(end)),'-len-',num2str(L0),'-eps',num2str(modelprops.epsilon),'-u',num2str(MV),'-a',num2str(modelprops.a)];
 
 L=L0*MV;
 %% IPE400
 h = (400)*10^(-3)*MV; %[m]
 b = 180e-3*MV; %[m]
 tw = 8.6e-3*MV; %[m]
 tf = 13.5e-3*MV; %[m]
 %% Load
 P = loadFactor*1e6*MV; %[N?]
 Emodul=2.1e+11/MV;
 %M = 0; % [N m?]
 
%% Finite Element Model
 
%  a=.1;

%  xcoords1a = ([linspace(-L,0,numofelm(1)+1)';linspace(0,L,numofelm(end)+1)']);
 xcoords1a = (linspace(-L,0,numofelm(1)+1)');
 xcoords1a(abs(xcoords1a)<1e-12) = 0;
%  xcoords1b=-xcoords1a;
%  xcoords2a=xcoords1a;
%  xcoords2b=-xcoords1a;
 ycoords1a =modelprops.a*xcoords1a;
%  ycoords1b=-ycoords1a;
%  ycoords2a=ycoords1a;
%  ycoords2b=-ycoords1a;

 %zcoords = 0*xcoords;

 NrNodi=length(xcoords1a);
 Nodes1a = [transpose(1:NrNodi), xcoords1a, ycoords1a];
 Nodes1b = [transpose(1:NrNodi)+NrNodi, -xcoords1a, -ycoords1a];
 Nodes2a = [transpose(1:NrNodi)+2*NrNodi, xcoords1a, -ycoords1a];
 Nodes2b = [transpose(1:NrNodi)+3*NrNodi, -xcoords1a, ycoords1a];
 NrEli=NrNodi-1;
 Elements1a = [transpose(1:NrEli),Nodes1a(1:end-1,1),Nodes1a(2:end,1)];
 Elements1b = [transpose(1:NrEli)+NrEli,Nodes1b(1:end-1,1),Nodes1b(2:end,1)];
 Elements2a = [transpose(1:NrEli)+2*NrEli,Nodes2a(1:end-1,1),Nodes2a(2:end,1)];
 Elements2b = [transpose(1:NrEli)+3*NrEli,Nodes2b(1:end-1,1),Nodes2b(2:end,1)];

 
 
 rp1a = Nodes1a(1,1);
 rp1b = Nodes1b(1,1);
 rp2a = Nodes2a(1,1);
 rp2b = Nodes2b(1,1);
 
 enableBeams=[-1 2];
 if numel(enableBeams)==4
  Nodes=[Nodes1a;Nodes1b;Nodes2a;Nodes2b];
  Elements=[Elements1a;Elements1b;Elements2a;Elements2b];
  rpall=[rp1a;rp1b;rp2a;rp2b];
 else
  Nodes=Nodes1a;
  Elements=Elements1a;
  rpall=rp1a;
  if any(enableBeams==1)
   Nodes=[Nodes;Nodes1b];
   Elements=[Elements;Elements1b];
   rpall=[rpall;rp1b];
  end
  if any(enableBeams==-2)
   Nodes=[Nodes;Nodes2a];
   Elements=[Elements;Elements2a];
   rpall=[rpall;rp2a];
  end
  if any(enableBeams==2)
   Nodes=[Nodes;Nodes2b];
   Elements=[Elements;Elements2b];
   rpall=[rpall;rp2b];
  end
  %   Nodes(:,1)=1:size(Nodes,1);
  rpallMat=rpall;
  for i=1:numel(rpall)
   rpallMat(i)=find(rpall(i)==Nodes(:,1));
  end
 end
 
 rpMiddle= Nodes(Nodes(:,2)==0,1);
 
 if strcmpi('B32',eltype(1:3)) || strcmpi('B22',eltype(1:3))
  Nodes2 = 0.5*(Nodes(1:end-1,2:end) + Nodes(2:end,2:end));
  Nodes2 = [Nodes(end,1) + ctranspose(1:size(Nodes2,1)),Nodes2];
  Nodes = [Nodes; Nodes2];
  Elements = [Elements(:,1),Elements(:,2),Nodes2(:,1),Elements(:,3)];%Nr,Anfang,mitte,Ende
 end

%    midnode1 = Nodes(Nodes(:,2)==0,1);
   rpmidnode1a=rpMiddle(1);
%    if numel(rpMiddle)>1
%     rpmidnode1b=rpMiddle(2);
%     if numel(rpMiddle)>2
%      rpmidnode2a=rpMiddle(3);
%      rpmidnode2b=rpMiddle(4);
%     end
%    end
%    Elements(Elements(:,2)==midnode1,2) = midnode2;
 
 %% Boundary conditions
 if strcmp(eltype,'B21') || strcmp(eltype,'B23') || strcmp(eltype,'B22')
  dofpNode=6;
 else
  error('MyProgram:Element','unknown Element')
 end
 %+1...Nx
 %+2...Ny
 %+6...Mz
 BCrpall=[1,2];
 BCrpMiddle=1;
%  BCmidnode2=BCrpMiddle;
%  BCrpRight=BCrpall;
 BC = ...
 [...
%   dofpNode*(rpall  - 1) + 1, 0;
%   dofpNode*(rpall  - 1) + 2, 0;
%   dofpNode*(rpMiddle- 1) + 2, 0;
%   dofpNode*(rpMiddle- 1) + 6, 0;
%   dofpNode*(rpRight - 1) + 2, 0;
  ];
 for i=1:numel(BCrpall)
  for j=1:numel(rpall)
   BC=[BC;dofpNode*(rpallMat(j) - 1) + BCrpall(i), 0]; %#ok<AGROW>
  end
 end
 for i=1:numel(BCrpMiddle)
  for j=1:numel(rpmidnode1a)
  BC=[BC;dofpNode*(rpmidnode1a(j) - 1) + BCrpMiddle(i), 0]; %#ok<AGROW>
  end
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
%   fprintf(u1,[num2str(transpose(rpall)),'\n']);
  fprintf(u1,'%d\n',rpall);
  fprintf(u1,'*Nset, nset=BCMitte, instance=Part-1-1\n');%for BC
%   fprintf(u1,[num2str(transpose(rpMiddle)),', ',num2str(transpose(midnode2)),'\n']);
  fprintf(u1,'%d\n',rpmidnode1a);
  fprintf(u1,'*Nset, nset=ALLMITTE, instance=Part-1-1\n');%for BC
  fprintf(u1,'%d\n',rpMiddle);
%   fprintf(u1,'*Nset, nset=MomemntMitte, instance=Part-1-1\n');
%   fprintf(u1,'%d, %d, %d, %d\n',rpMiddle);
  fprintf(u1,'*Nset, nset=FORCEMITTE, instance=Part-1-1\n');%for Force
%   fprintf(u1,[num2str(transpose(midnode2)),'\n']);
  fprintf(u1,'%d\n',rpmidnode1a);
%   fprintf(u1,'*Nset, nset=rightend, instance=Part-1-1\n');
%   fprintf(u1,'%d\n',rpall);
  
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+1,rpmidnode1a,rpmidnode1b]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+2,rpmidnode1a,rpmidnode2a]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+3,rpmidnode1a,rpmidnode2b]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');

%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+4,rpmidnode1b,rpmidnode2a]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+5,rpmidnode1b,rpmidnode2b]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');
%   fprintf(u1,'*Element,Type=CONN3D2, Elset=connector\n');
%   fprintf(u1,'%d, Part-1-1.%d, Part-1-1.%d\n',[size(Elements,1)+6,rpmidnode2a,rpmidnode2b]);
%   fprintf(u1,'*Connector Section, Elset=connector\n');
%   fprintf(u1,'JOIN\n');

fprintf(u1,'*Surface, type=NODE, name=ALLMITTE_CNS_, internal\n');
fprintf(u1,'ALLMITTE, 1.\n');
fprintf(u1,'** Constraint: Constraint-1\n');
fprintf(u1,'*Coupling, constraint name=Constraint-1, ref node=FORCEMITTE, surface=ALLMITTE_CNS_\n');
fprintf(u1,'*Kinematic\n');
fprintf(u1,'1, 1\n');
fprintf(u1,'2, 2\n');
  fprintf(u1,'*End Assembly\n');
  
  fprintf(u1,'*Material, name=Material-1\n');
  fprintf(u1,'*Elastic\n');
  fprintf(u1,[num2str(Emodul),', 0.3\n']);
  
  %% Boundary conditions
  %%
  
  fprintf(u1,'*Boundary\n');
%   %fprintf(u1,'leftend,  1, 1\n');
%   %fprintf(u1,'leftend,  2, 2\n');
%   fprintf(u1,'leftend,  3, 3\n');
%   fprintf(u1,'leftend,  4, 4\n');
%   %fprintf(u1,'rightend,  1, 1\n');
%   %fprintf(u1,'rightend,  2, 2\n');
%   fprintf(u1,'rightend,  3, 3\n');
%   fprintf(u1,'rightend,  4, 4\n');
 for i=1:numel(BCrpall)
  FoN=BCrpall(i);
  fprintf(u1,'leftend,  %d, %d\n',FoN,FoN);
 end
 for i=1:numel(BCrpMiddle)
  FoN=BCrpMiddle(i);
  fprintf(u1,'BCMitte,  %d, %d\n',FoN,FoN);
 end
%  for i=1:numel(BCrpRight)
%   FoN=BCrpRight(i);
%   fprintf(u1,'rightend,  %d, %d\n',FoN,FoN);
%  end
  
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
   fprintf(u3,'FORCEMITTE, 2, %f\n',-lambda(k)*P);
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
