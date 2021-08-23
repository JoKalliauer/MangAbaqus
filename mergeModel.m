function model=mergeModel(modelP1,modelM1)

model=modelP1;

model.lambda=[modelP1.lambda;NaN;-modelM1.lambda];
model.fulllambda=[modelP1.fulllambda;NaN;-modelM1.fulllambda];
[Zeilen , ~]=size(modelP1.fullEV);
model.fullEV=[modelP1.fullEV,NaN(Zeilen,1),modelM1.fullEV];
model.DetKtx=[modelP1.DetKtx;NaN;modelM1.DetKtx];
model.arclengthuJK=[modelP1.arclengthuJK;NaN;modelM1.arclengthuJK];
model.dxidl=[modelP1.dxidl;NaN;modelM1.dxidl];
model.lambda0=[modelP1.lambda0;NaN;-modelM1.lambda0];
model.load=[modelP1.load;NaN;-modelM1.load];
model.load0=[modelP1.load0;NaN;-modelM1.load0];
model.fullload0=[modelP1.fullload0;NaN;-modelM1.fullload0];
model.eigenvalues=[modelP1.eigenvalues;NaN*modelP1.eigenvalues{1};modelM1.eigenvalues];
NanEV=NaN*modelP1.eigenvectors{1};
if numel(NanEV)==0
 model.eigenvectors=[];
 model.eigvecDRH=[];
else
 model.eigenvectors=[modelP1.eigenvectors;NaN*modelP1.eigenvectors{1};modelM1.eigenvectors];
 model.eigvecDRH=[modelP1.eigvecDRH;NaN*modelP1.eigvecDRH{1};modelM1.eigvecDRH];
end
model.arclengths=[modelP1.arclengths;NaN*modelP1.arclengths{1};modelM1.arclengths];
model.arclengthurHM=[modelP1.arclengthurHM;NaN*modelP1.arclengthurHM{1};modelM1.arclengthurHM];
model.arclengthuHM=[modelP1.arclengthuHM;NaN*modelP1.arclengthuHM{1};modelM1.arclengthuHM];
model.uMaxJK=[modelP1.uMaxJK;NaN*modelP1.uMaxJK{1};modelM1.uMaxJK];
end

