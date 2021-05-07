function [R_DRHi] = RealEV(bungleEV,BC,R_DRHsize,nA)
 %inNr=inDOF(2)-inDOF(1)+1; % add number of internal (inverse) Nodes)
 R_DRHsizeIn=R_DRHsize(2:4);%+[0 inNr 0]; %DoFpNode x Nodes(inkl.Hyb) x NrEigs
 R_DRHi=NaN(R_DRHsizeIn);
 R_DRHij=NaN(R_DRHsizeIn(1:2));
%  ReducedHybridDofa=inDOF(3)+1:R_DRHsizeIn(1);
%  ReducedHybridDofb=inDOF(4)+1:R_DRHsizeIn(1);
%  nA=false(R_DRHsizeIn(1:2));
%  nA(ReducedHybridDofa,inDOF(1):2:inDOF(2))=true;
%  nA(ReducedHybridDofb,inDOF(1)+1:2:inDOF(2))=true;
 pos=1:prod(R_DRHsizeIn(1:2));
 Apos=pos(~nA);%Active Positions
 nBC=Apos;
 nBC(BC)=[];%not Boundary
 for j=1:size(bungleEV,2)
  bungleEVj=bungleEV(:,j);
  R_DRHij(BC)=0;
  R_DRHij(nA)=0;
  R_DRHij(nBC)=bungleEVj;
%   if ~any(isnan(R_DRHij(:))); assert(abs(norm(R_DRHij(:))-1)<.01,'norm not 1')end
  R_DRHi(:,:,j)=R_DRHij;
 end
end