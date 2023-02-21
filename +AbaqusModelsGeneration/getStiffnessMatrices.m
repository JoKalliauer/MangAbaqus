function [StifMatrices,num0,activeDofs,BC,inDOF,dofs] = getStiffnessMatrices(model,lambdareq,typeofanalysis)
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
 fname = [model.AbaqusRunsFolder,filename,'_STIF',num2str(num(i)),'.mtx'];
 %disp(fname);
 if usejava('jvm')
  %wb.Children.Title.Interpreter = 'none'; %https://de.mathworks.com/matlabcentral/answers/78895-how-do-i-make-interpreter-none-work-inside-the-waitbar-text#answer_298200
  waitbar(i/maxsteps,wb,fname,'interpreter','none');
 end
 mtxSparse = AbaqusModelsGeneration.translateStiffnessMtxFormatFromAbq(fname);
 dofs = max(mtxSparse(:,1:4)); % NrNodes(notHybrid), DegreesOfFreedomPerNode(inklWrapping), Nodes(notHybrid), DegreesOfFreedomPerNode(inklWrapping)

 %avoid internat degrees of freedom
 internalNodes = unique(mtxSparse(mtxSparse(:,1)<0,1));%list all negative Nodes
 internalNodesNewFirst=1+dofs(1);
 internalNodesNewLast=length(internalNodes)+dofs(1);
 internalNodesNew=ctranspose(internalNodesNewFirst:internalNodesNewLast);%internalNodesNew = ctranspose(1:length(internalNodes)) + dofs(1);% assign negative Nodes a pos Node-Value
 if i==1
  dofsinv = min(mtxSparse(:,1:4)); %#ok<NASGU>
  NrInNod=numel(internalNodes);
  if NrInNod>0
   inDOFpNList=NaN(NrInNod,1);
   for j=1:NrInNod
    inDOFpNList(j) = max(mtxSparse(mtxSparse(:,1)==internalNodes(j),2));
   end
   inDOFpNa = max(mtxSparse(mtxSparse(:,1)==internalNodes(1),2));
   assert(inDOFpNa==6);
   if NrInNod>1
    inDOFpNb = max(mtxSparse(mtxSparse(:,1)==internalNodes(2),2));
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
 end
 for j = 1:length(internalNodes)
  mtxSparse(mtxSparse(:,1)==internalNodes(j),1) = internalNodesNew(j);
  mtxSparse(mtxSparse(:,3)==internalNodes(j),3) = internalNodesNew(j);
 end
 
 %activeNodes = unique(mtxSparse(:,1));
 
 if i==1 % initial stiffness matrix
  StifMatrices{i,1} = 0;
  dofslist = zeros(dofs(1)*dofs(2),1);
  BC2 = BCMatlab;
 else
  StifMatrices{i,1} = lambda(i-1);
 end
 
 mtx = zeros(size(mtxSparse,1),3);
 %n = 0;
 for j = 1:size(mtxSparse,1)
  row = (mtxSparse(j,1)-1)*dofs(2) + mtxSparse(j,2);
  col = (mtxSparse(j,3)-1)*dofs(4) + mtxSparse(j,4);
  %             n = n+1;
  mtx(j,1) = row;
  mtx(j,2) = col;
  mtx(j,3) = mtxSparse(j,5);
  if (row==col)&&(i==1)
   dofslist(row) = row;
   if mtx(j,3)==1e+36
    BCAbaqus=[BCAbaqus;row]; %#ok<AGROW>
   end
  end
  %                n = n+1;
  %                mtx(n,1) = col;
  %                mtx(n,2) = row;
  %                mtx(n,3) = mtxSparse(j,5);
  %             end
 end
 
 BCsort=unique(BCMatlab(:,1));
 if i==1
  assert(dofs(2)==dofs(4),'insconsistent freedoms')
  if sum(strcmp(fieldnames(model), 'dofpNode')) == 0
   warning('MyProgram:notExist','model.dofpNode does not exist')
  else
   assert(dofs(2)==model.dofpNode,'dofpNode %d does not agree with Abaqus %d',model.dofpNode,dofs(2))
  end
  assert(dofs(1)==dofs(3),'insconsistent nr of nodes')
  assert(dofs(1)==size(model.Nodes,1),'nr of nodes do not agree with Abaqus')
 end
 assert(numel(BCsort)==numel(BCAbaqus),'number of BC does not agree')
 if all(BCsort==BCAbaqus)
  BC=BCMatlab;
 else
  warning('MyProgram:Problem','Check BC')
  error('MyProgram:Problem','Check BC')
  BC=BCMatlab; %#ok<UNRCH>
 end

 
 %mtx2 = sparse(mtx(:,1),mtx(:,2),mtx(:,3));
 %diagmtx2 = speye(size(mtx2,1)).*diag(mtx2);
 %mtx2 = mtx2+mtx2'-diagmtx2;
 
 if i==1
  dofslist(dofslist==0) = [];
  %             for k = 1:size(BC,1)
  %                dofslist(dofslist==BC(k,1)) = [];
  %             end
  dofslist = [dofslist,ctranspose(1:length(dofslist))]; %#ok<AGROW>
  dofslist2 = zeros(max(dofslist(:,1)),2);
  dofslist2(dofslist(:,1),:) = dofslist;
 end
 %unactiveDofs = find(dofslist2(:,1)==0);
 %numofDofs = size(Nodes,1)*dofs(2);
 activeDofs = find(dofslist2(:,1)~=0);
%  if length(unactiveDofs)+length(activeDofs) < numofDofs
%   %unactiveDofs = [unactiveDofs; ctranspose((numofDofs-dofs(2)+1):numofDofs)];
%  end
 for k = 1:length(dofslist)
  BC(BC2(:,1)==dofslist(k,1),1) = dofslist(k,2);
 end
 
 for k = 1:size(mtx,1)
  mtx(k,1) = dofslist2(mtx(k,1),2);
  mtx(k,2) = dofslist2(mtx(k,2),2);
 end
 
 mtx = sparse(mtx(:,1),mtx(:,2),mtx(:,3));
 diagmtx = speye(size(mtx,1)).*diag(mtx);
 mtx = mtx+mtx'-diagmtx;
 
 %         mtx = applyStiffBeam(mtx,rpLeft,leftnodes,Nodes);
 %         mtx = applyStiffBeam(mtx,rpRight,rightnodes,Nodes);
 
 mtx(BC(:,1),:) = 0;
 mtx(:,BC(:,1)) = 0;
 mtx(BC(:,1),BC(:,1)) = inf*speye(size(BC,1),size(BC,1)); %1e36*speye(size(BC,1),size(BC,1));
 
 StifMatrices{i,2} = mtx;
end
if usejava('jvm'); waitbar(1,wb,'getStiffnessMatrices finish'); close(wb);end
end