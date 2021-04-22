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
model.eigenvalues=[modelP1.eigenvalues;NaN*modelP1.eigenvalues{1};modelM1.eigenvalues];
model.eigenvectors=[modelP1.eigenvectors;NaN*modelP1.eigenvectors{1};modelM1.eigenvectors];
model.arclengths=[modelP1.arclengths;NaN*modelP1.arclengths{1};modelM1.arclengths];
end

