function [evec0,eval0,numofeigs] = solveCLEforMinEigNew(Kt,Ktprim,Eigres,Kt0_0,typeofanal,iter,model,modelprops,Ktprim0_0)

% if ~exist('model','var')
%  model.filename='';
% end

 numofeigs = modelprops.numofeigs;
 LM=Kt;%Linke Matrix fuers EW-Problem
%  if numel(LM)>4000000
%   fast=true;
%  else
%   fast=false;
%  end
 sortJK=1;%default is sort for the smallest abs-eigenvalue
 
 if nargin==0
  %evec0 = []; eval0 = [];
  error('MyProgram:Strange','no input')
 end

 if sum(strcmp(fieldnames(modelprops), 'sigma')) == 0
  modelprops.sigma=0;
 end
 

  DOFKt=size(Kt,1);
  if numofeigs>size(Kt,1)
   numofeigs = size(Kt,1);
  end
  if iter<1
   iter=1;%avoiding errors
   Kt=0*Kt;%that it returns NaNs;
  end
  switch typeofanal
   case 'I'
    B = speye(size(Kt,1),size(Kt,1));
    %EWgesucht=0;
    sortJK=0;
   case 'CLE'
    B = Ktprim;
    %EWgesucht=0;
   case 'I-CLE'
    B = -speye(size(Kt,1),size(Kt,1)) + Ktprim;
    %EWgesucht=0;
   case 'I+CLE'
    B = speye(size(Kt,1),size(Kt,1)) + Ktprim;
    %EWgesucht=0;
   case 'test'
    B = (Kt - Kt0_0);
    %EWgesucht=0;
   case 'test2'
    B = (Kt - Kt0_0);
    LM = Kt0_0;
    %EWgesucht=0;
   case 'I-CLE-test3'
    B = -speye(size(Kt,1),size(Kt,1))*sum(diag(Ktprim))/size(Ktprim,1) + Ktprim;
    %EWgesucht=0;
   case 'I-CLE-test4'
    B = speye(size(Kt,1),size(Kt,1))*sum(diag(Ktprim))/size(Ktprim,1) - Ktprim;
    %EWgesucht=0;
   case 'diag(CLE)-CLE'
    B = -diag(Ktprim)+ Ktprim;
    %EWgesucht=0;
   case 'KNL' % [ (Kts+Ktu) - EW * Kt0 ] 
    LM = Kt - Kt0_0;
    if sum(full(LM(:)))==0
     LM = Kt0_0;
    end
    B = speye(size(LM,1),size(LM,1));
    %EWgesucht=-1;
   case 'KNL2' % [ Kt - EW * Kt0 ]
    if iter>1
     B =  -Kt0_0;
    else
     evec0 = NaN(DOFKt,numofeigs);
     eval0=ones(numofeigs,1);
     return
     %B = NaN(size(Kt,1),size(Kt,1)); %B = eye(size(Kt,1),size(Kt,1));
    end
    %EWgesucht=0;
    sortJK=-1;%take the most negativ real value
   case 'KNL3' % [ Kt0 + EW * (Kts+Ktu) ]
    %B = Kt - Kt0_0;
    LM = Kt0_0;
    B=NaN(size(LM));
    %EWgesucht=1.1;
    RM=-(Kt - Kt0_0);
    if full(sum(RM(:)))==0
     warning('MyProgram:Empty','B is zero')
     RM = NaN(size(LM,1),size(LM,1));
    end
   case 'KNL4' % [ Kt0 - EW * (Kts+Ktu) ]
    RM = Kt - Kt0_0;
    LM = Kt0_0;
    if iter==1
     RM = NaN(size(LM));% LM = eye(size(LM));
    end
    if full(sum(RM(:)))==0
     warning('MyProgram:Empty','B is zero')
     RM = speye(size(LM,1),size(LM,1));
    end
    B=NaN(size(LM));
    %EWgesucht=-1;
   case 'Kg'
    eigenvalues = Eigres.eigenvalues{iter};
    eigenvectors = Eigres.eigenvectors;
    numofeigs = length(eigenvalues);
    
    keys = eigenvectors.keys;
    
    test = eigenvectors(keys{1});
    numofnodes = size(test{1},1);
    ndof = length(keys)*numofnodes;
    evec0 = NaN(ndof,numofeigs);
    for i = 1:length(keys)
     k1 = keys{i};
     dis1 = eigenvectors(k1);
     for j = 1:numofeigs
      evec0(i:length(keys):ndof,j) = dis1{j}(:,iter+1);
     end
    end
    eval0 = eigenvalues;
    numbers = ctranspose(1:numofeigs);
    posit = numbers(eval0>=0);
    negat = numbers(eval0<0);
    eval0 = [eval0(posit);eval0(negat)];
    evec0 = [evec0(:,posit),evec0(:,negat)];
    for ki = 1:size(evec0,2)
     evec0(:,ki) = evec0(:,ki)/norm(evec0(:,ki));
    end
    return
   case 'KsigmaKt0'
    error('MyProgram:NotImplemented','typeofanal: %s is unknown',typeofanal)
   case 'KNoLinearKt0'
    if iter>2
     %B =  -Kt0_0;
     RM = Kt0_0;
    else
     evec0 = NaN(DOFKt,numofeigs);
     eval0=ones(numofeigs,1);
     return
     %B = NaN(size(Kt,1),size(Kt,1)); %B = eye(size(Kt,1),size(Kt,1));
    end
    %EWgesucht=0;
    LM=Kt-model.fulllambda(iter)*Ktprim0_0;
    sortJK=-1;%take the most negativ real value
   otherwise
    warning('MyProgram:NotImplemented','typeofanal: %s is unknown',typeofanal)
    error('MyProgram:NotImplemented','typeofanal: %s is unknown',typeofanal)
  end %switch typeofanal
  
  if full(sum(LM(:)))==0
   warning('MyProgram:Empty','LM is zero')
   evec0 = NaN(DOFKt,numofeigs);
   eval0=ones(numofeigs,1);
   return
  end
  
  if exist('B','var')
   if full(sum(B(:)))==0
    warning('MyProgram:Empty','B is zero')
    B = NaN(size(LM,1),size(LM,1));
   end
  end
  
  if ~exist('RM','var')%~strcmpi(typeofanal,'KNL4') && ~strcmpi(typeofanal,'KNL3')
   RM=-B;
  end

  
  if any(any(~isnan(RM))) && any(any(Kt))
    %iter %#ok<NOPRT>
    %     if iter==10
    %      [row,col,v] = find(RM);
    %      dlmwrite(strcat('Output/',model.filename,'B.txt'),[col,row,v], 'delimiter', '\t')
    %      [row,col,v] = find(LM);
    %      dlmwrite(strcat('Output/',model.filename,'Kt.txt'),[col,row,v], 'delimiter', '\t')
    %     end
    %     if EWgesucht==0
    %      if 1==0
    % %       %check if 0 and 'smallestreal' are the same
    % %       [evec0,eval0] = eigs(LM,RM,numofeigs,0);
    % %       [evecsr,evalsr] = eigs(LM,RM,numofeigs,'smallestreal');
    % %       diageval0 = diag(eval0);
    % %       diagevalsr = diag(evalsr);
    % %       if 1==1
    % %        keep=[true;(diageval0(2:end)==diagevalsr(2:end))];
    % %       else
    % %        keep=true(numofeigs,1);
    % %       end
    % %       if ~all(all(evec0(:,keep)==evecsr(:,keep)))
    % %        %firstProblem=find(~(diageval0==diagevalsr));
    % %        %warning('Myprogram:SRnot0','smallestreal differ from 0, first eigenvaluedifference is %i',firstProblem(1));
    % %       end
    % %       evec=evec0(:,keep);
    % %       eval=eval0(:,keep);
    % %       numofeigs=numofeigs-sum(~keep);
    %      end
    %     else
    %     end
    assert(any(any(LM~=RM)),'LM and RM must be different')
%     if fast==true
%      [evec,eval] = eigs(LM,RM,numofeigs,modelprops.sigma);
%     else
     [evec,eval] = eigs(LM,RM,numofeigs,modelprops.sigma,'Tolerance',eps(1e-290),'MaxIterations',3000);
%     end
    diageval = diag(eval);
    diffs=diageval(2:end)-diageval(1:end-1);
    relevantAbs=(abs(diffs)<2.923e-15);%1.0e-06
    relevantRel=(abs(diffs./diageval(1:end-1))<3.962e-10);%1.0e-06
    relevant=logical(relevantAbs+relevantRel); %sobald min einer 1 ist ist relevant ungleich 0
    if any(relevant)
     %warning('MyProgram:Precission','close eigenvalues, setting eigenvectors to NaN')
     %diageval([relevant;false])=NaN;
     %diageval([false;relevant])=NaN;
     evec(:,[relevant;false])=NaN;
     evec(:,[false;relevant])=NaN;
    end
    if ~isreal(eval)
     if iter<6 || any(isnan(diageval))
      evec0 = NaN(DOFKt,numofeigs);
      eval0=NaN(numofeigs,1);
      return
     end
     relations=abs(imag(diageval(:)))./(real(diageval(:)));
     relation=max(relations(1));
     %warning('MyProgram:NaN','komplexer Anteil des Eigenwertes ist %s',relation)
     %check of first eigenvalue
     if relation>0.0466%0.00481
      warning('MyProgram:NaN','komplexer Anteil des 1.Eigenwertes ist %s',relation);
      if relation>0.686
       error('MyProgram:complex','komplex Part too large')
       %evec0 = NaN(DOFKt,numofeigs);
       %eval0=NaN(numofeigs,1);
       %return
      end
     end
     
%      %check of higher eigenvalues
%      komplex = (relations>0.0466);
%      if any(komplex)%0.00481
%       warning('MyProgram:NaN','komplexer Anteil eines Eigenwertes ist %s',max(relations));
%      end
     
    end
    
    %         [evec,eval] = eig(full(Kt),-full(Ktprim));
  else %B is NaN
   warning('MyProgram:NaN','B is NaN or Kt=0');
   evec=NaN(size(RM,1),numofeigs);
   eval=NaN(numofeigs,numofeigs);
   diageval = diag(eval);
  end
  
 
  if sortJK==0
   % do not change order
   s=1:numofeigs;
  elseif sortJK==1
   %take the value closest to zero
   [~,s] = sort(abs(diageval));
  elseif sortJK==-1
   %take the most negativ real value
   [~,s] = sort(real(diageval));
  else
   error('MyProgram:Undefinded','modelprops.sort unknown')
  end
  evec0 = evec(:,s(1:numofeigs));
  evec0(abs(evec0)<8.5e-14)=0;
  eval0 = diageval(s(1:numofeigs));
  if ~isreal(eval0)
   warning('MyProgram:Komplex','eval0 is complex')
  end
  
  if modelprops.allowComplex==false
   relations=abs(imag(eval0(:)))./(real(eval0(:)));
   %check of higher eigenvalues
   komplex = (relations>0.0005);
   if any(komplex)%0.00481
    %warning('MyProgram:NaN','komplexer Anteil eines Eigenwertes ist %s',max(relations));
    evec0(:,komplex)=NaN;
    eval0(komplex,1)=NaN;
    if sortJK==0
     % do not change order
     s=1:numofeigs;
    elseif sortJK==1
     %take the value closest to zero
     [~,s] = sort(abs(eval0));
    elseif sortJK==-1
     %take the most negativ real value
     [~,s] = sort(real(eval0));
    else
     error('MyProgram:Undefinded','modelprops.sort unknown')
    end
    %numofeigs=max(sum(~isnan(eval0)),1);
    evec0 = evec0(:,s(1:numofeigs));
    eval0 = eval0(s(1:numofeigs));
   end
  
   for i=1:numofeigs
    if isreal(eval0(i))
     evec0(:,i)=evec0(:,i)/norm(evec0(:,i));
    else
     evec0(:,i)=NaN(size(evec0(:,i)));
    end
   end
   
   if abs(norm(evec0(:,1))-1)>eps(1)
    assert(abs(norm(evec0(:,1))-1)<=eps(1),'norm nicht eins');
   end
   
   numofeigs=find(~isnan(eval0), 1, 'last' );%currently no use of this line
  end


 for i = 1:size(evec0,2)
  evec0(:,i) = evec0(:,i)/norm(evec0(:,i));
 end
 %assert(isreal(eval0),'eval0 is complex')
end
