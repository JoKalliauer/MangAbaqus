function model=mergeModel(modelP1,modelM1)
% this function merges the solution for \lambda>0 with the solution of \lambda<0

%% Input
% modelP1 ... results with modelprops.loadfactor=+1;
% modelM1 ... results with modelprops.loadfactor=-1;

%% Output
% model has the results of modelP1 and modelM1 merged into one model.

%% Recent changes
% 2023-04-13 JK: if main.whichEV='skip' there might not be an eigenvector

model=modelP1;

model.lambda=[modelP1.lambda;NaN;-modelM1.lambda];
model.fulllambda=[modelP1.fulllambda;NaN;-modelM1.fulllambda];
[Zeilen , ~]=size(modelP1.fullEV);
model.fullEV=[modelP1.fullEV,NaN(Zeilen,1),modelM1.fullEV];
model.DetKtx=[modelP1.DetKtx;NaN;modelM1.DetKtx];
model.arclengthuJK=[modelP1.arclengthuJK;NaN;modelM1.arclengthuJK];
model.dxidl=[modelP1.dxidl;NaN;modelM1.dxidl];
model.lambdainput=[modelP1.lambdainput;NaN;-modelM1.lambdainput];
model.load=[modelP1.load;NaN;-modelM1.load];
model.load0=[modelP1.load0;NaN;-modelM1.load0];
model.fullload0=[modelP1.fullload0;NaN;-modelM1.fullload0];
model.eigenvalues=[modelP1.eigenvalues;NaN*modelP1.eigenvalues{1};modelM1.eigenvalues];
if sum(strcmp(fieldnames(modelP1), 'eigenvectors')) ~= 0 %skip it main.whichEV='skip'
 NanEV=NaN*modelP1.eigenvectors{1};
 if numel(NanEV)==0
  model.eigenvectors=[];
  model.eigvecDRH=[];
 else
  model.eigenvectors=[modelP1.eigenvectors;NanEV;modelM1.eigenvectors];
  model.eigvecDRH=[modelP1.eigvecDRH;NaN*modelP1.eigvecDRH{1};modelM1.eigvecDRH];
 end
end
model.arclengths=[modelP1.arclengths;NaN*modelP1.arclengths{1};modelM1.arclengths];
model.arclengthurHM=[modelP1.arclengthurHM;NaN*modelP1.arclengthurHM{1};modelM1.arclengthurHM];
model.arclengthuHM=[modelP1.arclengthuHM;NaN*modelP1.arclengthuHM{1};modelM1.arclengthuHM];
model.uMaxJK=[modelP1.uMaxJK;NaN*modelP1.uMaxJK{1};modelM1.uMaxJK];
end

