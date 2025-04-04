function [StifMatrices] = getStiffnessMatrices3(model,lambdareq,typeofanalysis)


if usejava('jvm'); wb=waitbar(0,'getStiffnessMatrices','title','getStiffnessMatrices3');end
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
cd ~/ownCloud/Post/MangAbaqus/ %cd ..

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
% if steps~=length(num)
%  warning('MyProgram:Inputchange','Abaqus calculated more Loadsteps than requested, try using modelprops.forceAbaqus=true')
%  num0 = num(1:steps);
% else
%  num0 = num;
% end
stepsAbaqus=length(num);
%num0Abaqus=num;
maxsteps=max(steps,stepsAbaqus);
BCAbaqus=[];
for i = 1:steps
 if usejava('jvm'); waitbar(i/maxsteps,wb,'getStiffnessMatrices ...');end
 fname = [model.AbaqusRunsFolder,filename,'_STIF',num2str(num(i)),'.mtx'];
 %disp(fname);
 if usejava('jvm')
  wb.Children.Title.Interpreter = 'none'; %https://de.mathworks.com/matlabcentral/answers/78895-how-do-i-make-interpreter-none-work-inside-the-waitbar-text#answer_298200
  waitbar(i/maxsteps,wb,fname,'interpreter','none');
 end
 mtxSparse = AbaqusModelsGeneration.translateStiffnessMtxFormatFromAbq(fname);
 model.dofpNode=2;
 mtxSparse(mtxSparse(:,2)>3,:)=[];
 mtxSparse(mtxSparse(:,4)>3,:)=[];
 mtxSparse(mtxSparse(:,2)<1,:)=[];
 mtxSparse(mtxSparse(:,4)<1,:)=[];
 dofs = max(mtxSparse(:,1:4));
 %dofsinv = min(mtxSparse(:,1:4)); %#ok<NASGU>
 %avoid internat degrees of freedom
 internalNodes = unique(mtxSparse(mtxSparse(:,1)<0,1));
 internalNodesNew = ctranspose(1:length(internalNodes)) + dofs(1);
 for j = 1:length(internalNodes)
  mtxSparse(mtxSparse(:,1)==internalNodes(j),1) = internalNodesNew(j);
  mtxSparse(mtxSparse(:,3)==internalNodes(j),3) = internalNodesNew(j);
 end
 
 
 if i==1 % initial stiffness matrix
  StifMatrices{i,1} = 0;
  dofslist = zeros(dofs(1)*dofs(2),1);
  BC2 = BCMatlab;
 else
  StifMatrices{i,1} = lambda(i-1);
 end
 
 mtx = zeros(size(mtxSparse,1),3);
 for j = 1:size(mtxSparse,1)
  row = (mtxSparse(j,1)-1)*dofs(2) + mtxSparse(j,2);
  col = (mtxSparse(j,3)-1)*dofs(4) + mtxSparse(j,4);
  mtx(j,1) = row;
  mtx(j,2) = col;
  mtx(j,3) = mtxSparse(j,5);
  if (row==col)&&(i==1)
   dofslist(row) = row;
   if mtx(j,3)==1e+36
    BCAbaqus=[BCAbaqus;row]; %#ok<AGROW>
   end
  end

 end
 
 BCsort=unique(BCMatlab(:,1));
 if i==1
  assert(dofs(2)==dofs(4),'insconsistent freedoms')
  assert(dofs(2)==model.dofpNode,'dofpNode %d does not agree with Abaqus %d',model.dofpNode,dofs(2))
  assert(dofs(1)==dofs(3),'insconsistent nr of nodes')
  assert(dofs(1)==size(model.Nodes,1),'nr of nodes do not agree with Abaqus')
 end
%  if all(BCsort==BCAbaqus)
%   BC=BCMatlab;
%  else
%   warning('MyProgram:Problem','Check BC')
%   error('MyProgram:Problem','Check BC')
%   BC=BCMatlab; %#ok<UNRCH>
%  end

 

 
 if i==1
  dofslist(dofslist==0) = [];

  dofslist = [dofslist,ctranspose(1:length(dofslist))]; %#ok<AGROW>
  dofslist2 = zeros(max(dofslist(:,1)),2);
  dofslist2(dofslist(:,1),:) = dofslist;
 end

%  activeDofs = find(dofslist2(:,1)~=0);

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
 

 
 mtx(BC(:,1),:) = 0;
 mtx(:,BC(:,1)) = 0;
 mtx(BC(:,1),BC(:,1)) = 1e36*speye(size(BC,1),size(BC,1));
 
 StifMatrices{i,2} = mtx;
end
if usejava('jvm'); waitbar(1,wb,'getStiffnessMatrices finish'); close(wb);end
end