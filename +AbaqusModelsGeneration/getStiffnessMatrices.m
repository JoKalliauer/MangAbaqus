function [StifMatrices,num0,activeDofs,BC,inDOF,dofs] = getStiffnessMatrices(model,lambdareq,typeofanalysis,elementtype)
%% get the stiffnessMatrices from Abaqus
%university:TU Wien
%author:Michał Malendowski (©2019-2020), Johannes Kalliauer(©2020-2023)


%% Input
% model ... struture containing data like filename or Abaqus-folder
% lambdareq
% typeofanalysis .. if ~strcmp(typeofanalysis,'Kg'), than some file don't need to exist


%% Output
% StifMatrices ... thats the stiffnessmatrix read out from Abauqusfiles
% num0 .. numbers of the Stiff-files
% activeDofs ... dofs that are not resticted like BoundaryConditions
% BC ... Location and value of Boundary-condition



%% Recent Changes
%2023-02-21 JK: added explantations for dofs
%2023-03-16 JK: assert(inDOFpNa==6); leaded to an error for TL_arch3D-B33-20-f1-eps0.005-u1-E210000000000-KNL2-1
%2023-03-31 JK: added comments for Stiffnesmatrix, deleted outcommented lines
%2023-03-31 JK: removed that the rows of the boundary conditions are set to zero

%% input check
%EulerBernulli=false;
if ~exist('elementtype','var')
 elementtype='unknown';
 isB32=true; %assuming it is one of the B32 elements
 elNr=NaN;
elseif strcmp(elementtype,'B33') || strcmp(elementtype,'B33H') || strcmp(elementtype,'B31H')
 isB32=false;
 elNr=33;
elseif strcmp(elementtype(2:3),'21') || strcmp(elementtype,'B21H') 
 elNr=21;
else
 isB32=true;
 elNr=32;
end

%% Code

if usejava('jvm'); wb=waitbar(0,'getStiffnessMatrices','name','getStiffnessMatrices','WindowState','minimized');end
%Loads = [];

filename = model.filename;
lambda = model.lambda;
BCMatlab = model.BCMatlab;
%Nodes = model.Nodes;

if ~exist(model.AbaqusRunsFolder, 'dir')
 if isunix
  warning('MyProgram:OS','AbaqusRunsFolder does not exist')
  mkdir(model.AbaqusRunsFolder);
 end
 if ispc
  warning('MyProgram:OS','You are using Windows and AbaqusRunsFolder does not exist, therfore skipping')
  return
 end
end
if ~exist([model.AbaqusRunsFolder,filename,'_STIF9.mtx'],'file')
 %warning('MyProgramm:Missing','_STIF*.mtx missing')
 if usejava('jvm'); waitbar(1,wb,'getStiffnessMatrices error');end
 if ~exist([model.AbaqusRunsFolder,filename,'_STIF7.mtx'],'file') && ~strcmp(typeofanalysis,'Kg')
  if usejava('jvm'); close(wb);end
  model.AbaqusRunsFolder
  error('MyProgramm:Missing','_STIF*.mtx missing in %s , try rerunning forceAbaqus=true',model.AbaqusRunsFolder)
  %return
 else
  warning('MyProgram:Abaqus','only few stif existing, maybe abaqus failed?')
 end
end
oldpwd=pwd;
cd(model.AbaqusRunsFolder) %cd /home/jkalliau/Abaqus/Post/MangAbaqus/AbaqusRuns %cd AbaqusRuns
lenFN = length(filename);
StifFilesMatlabVersion = ls([filename,'_STIF*.mtx']);
StifFilesMatlabVersionTranspose=StifFilesMatlabVersion.';
StifFiles = convertCharsToStrings(StifFilesMatlabVersionTranspose);
beg_ = strfind(StifFiles,[filename,'_STIF']);
end_ = strfind(StifFiles,'.mtx');
num = zeros(length(beg_),1);
for i = 1:length(beg_)
 num(i) = str2double(StifFilesMatlabVersionTranspose((beg_(i)+lenFN+5):(end_(i)-1)));
end
cd(oldpwd);%cd ~/ownCloud/Post/MangAbaqus/ %cd ..

num = sort(num);

%     lambda = lambda(1:(length(num)-1));
if ~exist('lambdareq','var')
 StifMatrices = cell(length(num),2);
elseif isempty(lambdareq)
 StifMatrices = cell(length(num),2);
else
 StifMatrices = cell(1,2);
 num = num([0;lambda]==lambdareq);
end
steps=min(length(num),length(lambda)+1);
if steps~=length(num)
 assert(length(num)>=steps,'Maybe wrong')
 if steps+3<length(num)
  warning('MyProgram:Inputchange','Abaqus calculated more Loadsteps than requested, try using modelprops.forceAbaqus=true')
 end
 num0 = num(1:steps);
else
 num0 = num;
end
stepsAbaqus=length(num);
%num0Abaqus=num;
maxsteps=max(steps,stepsAbaqus);
BCAbaqus=[];
for i = 1:steps
 if usejava('jvm'); waitbar(i/maxsteps,wb,'getStiffnessMatrices ...');end
 fname_i = [model.AbaqusRunsFolder,filename,'_STIF',num2str(num(i)),'.mtx'];
 if usejava('jvm')
  waitbar(i/maxsteps,wb,fname_i,'interpreter','none');
 end
 mtxSparse_i = AbaqusModelsGeneration.translateStiffnessMtxFormatFromAbq(fname_i);
 if i==1
  dofs = max(mtxSparse_i(:,1:4)); % NrNodes(notHybrid), DegreesOfFreedomPerNode(inklWrapping), Nodes(notHybrid), DegreesOfFreedomPerNode(inklWrapping)
  
  %avoid internat degrees of freedom
  internalNodes = unique(mtxSparse_i(mtxSparse_i(:,1)<0,1));%list all negative Nodes
  internalNodesNewFirst=1+dofs(1);
  internalNodesNewLast=length(internalNodes)+dofs(1);
  internalNodesNew=ctranspose(internalNodesNewFirst:internalNodesNewLast);%internalNodesNew = ctranspose(1:length(internalNodes)) + dofs(1);% assign negative Nodes a pos Node-Value
  
  dofsinv = min(mtxSparse_i(:,1:4)); %#ok<NASGU>
  NrInNod=numel(internalNodes);
  if NrInNod>0
   inDOFpNList=NaN(NrInNod,1);
   for j=1:NrInNod
    inDOFpNList(j) = max(mtxSparse_i(mtxSparse_i(:,1)==internalNodes(j),2));
   end
   inDOFpNa = max(mtxSparse_i(mtxSparse_i(:,1)==internalNodes(1),2));
   if elNr==32
    assert(inDOFpNa==6);
   elseif elNr==33
    warning('MyPrgm:NoCheckImplemented','no check for B33 implemented')
    if strcmp(elementtype,'B33H') || strcmp(elementtype,'B31H')
     assert(inDOFpNa==3);%dont know why it is three
    elseif strcmp(elementtype,'B33')
     assert(inDOFpNa==2);%dont know why it is two
    else
     error('MyPrgm:Unknown','not implemted')
    end
   elseif elNr==21
    assert(inDOFpNa==2);%dont know why it is two
   else
    error('MyPrgm:Missing','check not implemented')
   end
   if NrInNod>1 && isB32
    inDOFpNb = max(mtxSparse_i(mtxSparse_i(:,1)==internalNodes(2),2));
    assert(inDOFpNb==6);
   else
    warning('MyPrgm:DontKnowIfItIsStrange','please remove this warning it might not be important')
    inDOFpNb=NaN;
   end
   inDOFpN=[inDOFpNa inDOFpNb];
  else
   inDOFpN=[0 0];
  end
  inDOF = [internalNodesNewFirst internalNodesNewLast inDOFpN];
 end %if i==1
 for j = 1:length(internalNodes)
  mtxSparse_i(mtxSparse_i(:,1)==internalNodes(j),1) = internalNodesNew(j);
  mtxSparse_i(mtxSparse_i(:,3)==internalNodes(j),3) = internalNodesNew(j);
 end
 
 %activeNodes = unique(mtxSparse(:,1));
 
 if i==1 % initial stiffness matrix
  StifMatrices{i,1} = 0;
  dofslist = zeros(dofs(1)*dofs(2),1);
  BC2 = BCMatlab;
 else
  StifMatrices{i,1} = lambda(i-1);
 end
 
 mtx_i = zeros(size(mtxSparse_i,1),3);
 for j = 1:size(mtxSparse_i,1) %loop over all entries in Sparsematrix
  row = (mtxSparse_i(j,1)-1)*dofs(2) + mtxSparse_i(j,2);
  col = (mtxSparse_i(j,3)-1)*dofs(4) + mtxSparse_i(j,4);
  mtx_i(j,1) = row;
  mtx_i(j,2) = col;
  mtx_i(j,3) = mtxSparse_i(j,5);
  if (row==col)&&(i==1)
   dofslist(row) = row;
   if mtx_i(j,3)==1e+36
    BCAbaqus=[BCAbaqus;row]; %#ok<AGROW>
   end
  end
 end %for j = 1:size(mtxSparse,1)
 
 if i==1
  BCsort=unique(BCMatlab(:,1));
  assert(dofs(2)==dofs(4),'insconsistent freedoms')
  if sum(strcmp(fieldnames(model), 'dofpNode')) == 0
   warning('MyProgram:notExist','model.dofpNode does not exist')
  else
   assert(dofs(2)==model.dofpNode,'dofpNode %d does not agree with Abaqus %d',model.dofpNode,dofs(2))
  end
  assert(dofs(1)==dofs(3),'insconsistent nr of nodes')
  assert(dofs(1)==size(model.Nodes,1),'nr of nodes do not agree with Abaqus')
  assert(numel(BCsort)==numel(BCAbaqus),'number of BC does not agree')  
  if all(BCsort==BCAbaqus)
   BC=BCMatlab;
  else
   warning('MyProgram:Problem','Check BC')
   error('MyProgram:Problem','Check BC')
   BC=BCMatlab; %#ok<UNRCH>
  end
  dofslist(dofslist==0) = [];
  dofslist = [dofslist,ctranspose(1:length(dofslist))]; %#ok<AGROW>
  dofslist2 = zeros(max(dofslist(:,1)),2);
  dofslist2(dofslist(:,1),:) = dofslist;
  activeDofs = find(dofslist2(:,1)~=0);
  for k = 1:length(dofslist)
   BC(BC2(:,1)==dofslist(k,1),1) = dofslist(k,2);
  end
 end % if i==1
 
 for k = 1:size(mtx_i,1)
  mtx_i(k,1) = dofslist2(mtx_i(k,1),2);
  mtx_i(k,2) = dofslist2(mtx_i(k,2),2);
 end
 
 mtx_i = sparse(mtx_i(:,1),mtx_i(:,2),mtx_i(:,3)); %lower trianular matrix
 diagmtx = speye(size(mtx_i,1)).*diag(mtx_i); %temporary value for mirrowing the matrix
 mtx_i = mtx_i+mtx_i'-diagmtx; % adding the values for the upper trianular
 
 %mtx_i(BC(:,1),:) = 0;
 %mtx_i(:,BC(:,1)) = 0;
 mtx_i(BC(:,1),BC(:,1)) = inf*speye(size(BC,1),size(BC,1)); %1e36*speye(size(BC,1),size(BC,1));
 
 StifMatrices{i,2} = mtx_i;
end %for i = 1:steps
if usejava('jvm'); waitbar(1,wb,'getStiffnessMatrices finish'); close(wb);end
end