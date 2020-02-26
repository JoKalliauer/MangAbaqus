function [StifMatrices,num0,activeNodes,activeDofs,unactiveDofs,BC,Loads] = getStiffnessMatrices(model,lambdareq)

Loads = [];

filename = model.filename;
lambda = model.lambda;
BC = model.BC;
Nodes = model.Nodes;
        
    cd AbaqusRuns
      lenFN = length(filename);
      StifFiles = ls([filename,'_STIF*.mtx']);
      beg_ = strfind(StifFiles,[filename,'_STIF']);
      end_ = strfind(StifFiles,'.mtx');
      num = zeros(length(beg_),1);
      for i = 1:length(beg_)
         num(i) = str2double(StifFiles((beg_(i)+lenFN+5):(end_(i)-1)));
      end
    cd ..
      
    num = sort(num);
    num0 = num;
%     lambda = lambda(1:(length(num)-1));
    if nargin<2
        StifMatrices = cell(length(num),2);
    elseif isempty(lambdareq)
        StifMatrices = cell(length(num),2);
    else
        StifMatrices = cell(1,2);
        num = num([0;lambda]==lambdareq);
    end
    for i = 1:length(num)
        fname = ['AbaqusRuns/',filename,'_STIF',num2str(num(i)),'.mtx'];
        disp(fname);
        mtxSparse = AbaqusModelsGeneration.translateStiffnessMtxFormatFromAbq(fname);
        dofs = max(mtxSparse(:,1:4));
        %avoid internat degrees of freedom
           internalNodes = unique(mtxSparse(mtxSparse(:,1)<0,1));
           internalNodesNew = [1:length(internalNodes)]' + dofs(1);
           for j = 1:length(internalNodes)
               mtxSparse(mtxSparse(:,1)==internalNodes(j),1) = internalNodesNew(j);
               mtxSparse(mtxSparse(:,3)==internalNodes(j),3) = internalNodesNew(j);
           end
        
        activeNodes = unique(mtxSparse(:,1));
         
        if i==1 % initial stiffness matrix
            StifMatrices{i,1} = 0;
            dofslist = zeros(dofs(1)*dofs(2),1);
            BC2 = BC;
        else
            StifMatrices{i,1} = lambda(i-1);
        end
        
        mtx = zeros(size(mtxSparse,1),3);
        n = 0;
        for j = 1:size(mtxSparse,1)
            row = (mtxSparse(j,1)-1)*dofs(2) + mtxSparse(j,2);
            col = (mtxSparse(j,3)-1)*dofs(4) + mtxSparse(j,4);
%             n = n+1;
            mtx(j,1) = row;
            mtx(j,2) = col;
            mtx(j,3) = mtxSparse(j,5);
            if (row==col)&&(i==1)
                dofslist(row) = row;
            end
%                n = n+1; 
%                mtx(n,1) = col;
%                mtx(n,2) = row;
%                mtx(n,3) = mtxSparse(j,5);
%             end
        end

        mtx2 = sparse(mtx(:,1),mtx(:,2),mtx(:,3));
        diagmtx2 = speye(size(mtx2,1)).*diag(mtx2);
        mtx2 = mtx2+mtx2'-diagmtx2;
        
        if i==1
            dofslist(dofslist==0) = [];
%             for k = 1:size(BC,1)
%                dofslist(dofslist==BC(k,1)) = []; 
%             end
            dofslist = [dofslist,[1:length(dofslist)]'];
            dofslist2 = zeros(max(dofslist(:,1)),2);
            dofslist2(dofslist(:,1),:) = dofslist;
        end
        unactiveDofs = find(dofslist2(:,1)==0);
        numofDofs = size(Nodes,1)*dofs(2);
        activeDofs = find(dofslist2(:,1)~=0);
        if length(unactiveDofs)+length(activeDofs) < numofDofs
            unactiveDofs = [unactiveDofs; [(numofDofs-dofs(2)+1):numofDofs]'];
        end
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
        mtx(BC(:,1),BC(:,1)) = 1e36*speye(size(BC,1),size(BC,1));
        
        StifMatrices{i,2} = mtx;
    end

end