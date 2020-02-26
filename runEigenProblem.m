function model = runEigenProblem(modelprops)
%modelname,lambda,epsil,numofelm,typeofanal,additionalParameters
if nargin<1
    testcase = 'TL_arch';
    numofelm = 20;
    eltype = 'B21';
    lambda = [.025,.075,.125,.175,.225,.2750,.325,.375,.425,.475,.525,.575,.625,.675,.725];
    epsil = 0.005;
    typeofanal = 'K0';
    loadFactor = 1.0;
    len = [];
else
    testcase = modelprops.testcase;
    numofelm = modelprops.numofelm;
    eltype =   modelprops.elementtype;
    lambda =   modelprops.lambda;
    epsil =    modelprops.epsilon;
    typeofanal = modelprops.typeofanalysis;
    loadFactor = modelprops.loadfactor;
    len = [];
end

if isempty(len)
    len = 0;
end

lambda = reshape(lambda,1,length(lambda));
    lambda0 = [0,lambda];
    lambda11 = lambda0 + epsil;
    lambda12 = lambda0 + 2*epsil;
    lambda13 = lambda0 + 3*epsil;
    lambda14 = lambda0 + 4*epsil;
    lambda15 = []; %lambda0 + 5*epsil;
    lambda21 = lambda - epsil;
    lambda22 = lambda - 2*epsil;
    lambda23 = lambda - 3*epsil;
    lambda24 = lambda - 4*epsil;
    lambda25 = []; %lambda - 5*epsil;
    lambda = [lambda0,lambda11,lambda21,lambda12,lambda22,lambda13,lambda23,lambda14,lambda24,lambda15,lambda25]';
    lambda = sort(unique(round(lambda*100000))/100000);
    lambda0 = sort(unique(round(lambda0*100000))/100000);

modelprops.lambda = lambda;
model = selectModel(modelprops);

filename = model.filename;
if exist(['AnalysisResults/',filename,'-',typeofanal,'.mat'], 'file') == 2
    load(['AnalysisResults/',filename,'-',typeofanal,'.mat'],'model');
    return
end

AbaqusModelsGeneration.runAbaqus(model.filename);
% Kg = AbaqusModelsGeneration.getKgmatrix(model);

[Kts,num,activeNodes,activeDofs,unactiveDofs,BC,Loads]  = AbaqusModelsGeneration.getStiffnessMatrices(model);
model.BC = BC;

       for j = 1:size(BC,1)
          activeDofs(activeDofs==BC(j,1)) = []; unactiveDofs = [unactiveDofs; BC(j,1)]; 
       end
       unactiveDofs = sort(unique(unactiveDofs));

[ELres, Nres, EigRes] = AbaqusModelsGeneration.getHistoryOutputFromDatFile(['AbaqusRuns/',model.filename,'.dat']);
% [membrane, nonmembrane] = AbaqusModelsGeneration.GetEnergies(ELres,model.Nodes,model.Elements);

Displ = NodalResults2Displ(Nres);
Kg = EigRes;

eigval = cell(length(lambda0),1);
eigvec = cell(length(lambda0),1);
displacements = cell(length(lambda0),1);
arclengths = cell(length(lambda0),1);
StiffMtxs = cell(length(lambda0),3);

matches = [];
n = 0;
m = 0;
for i = 1:(size(num,1)-1)
    if ~isempty(find(lambda0==lambda(i)))
        n = n+1;
        matches(n) = i;
    end
end
if matches(1)~=1
    matches = [1,matches];
end

Energy = zeros(length(lambda0),3);
Energy(:,1) = lambda0;
for i = 1:length(matches)
    disp('Lambda:');
    disp(lambda(matches(i)));
    disp('-------------');
%     Energy(i,2) = membrane(matches(i));
%     Energy(i,3) = nonmembrane(matches(i));
    if i==1
       [~,~,numofeigs] = solveCLEforMinEigNew();
       % 0.
       Kt0 = Kts{matches(i),2};
       rr = zeros(size(Kt0,1),1);
       ra = [1:size(Kt0,1)]';
%        ru = [];
%        ru = BC(:,1);
       ru = diag(Kt0==1e36);
       ra(ru) = [];
       RR0 = zeros(size(Kt0,1),numofeigs);
       Kt0(ru,:) = []; Kt0(:,ru) = [];
       Kt0_0 = Kt0;
%        Kg(ru,:) = []; Kg(:,ru) = [];
       % 1.
       dksi11 = sqrt(Displ{matches(i)}'*Displ{matches(i)});
       Kt11 = Kts{matches(i)+1,2}; Kt11(ru,:) = []; Kt11(:,ru) = [];
       % 2.
       dksi12 = sqrt((Displ{matches(i)+1}-Displ{matches(i)})'*(Displ{matches(i)+1}-Displ{matches(i)}));
       Kt12 = Kts{matches(i)+2,2}; Kt12(ru,:) = []; Kt12(:,ru) = [];
       % 3.
       dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));
       Kt13 = Kts{matches(i)+3,2};  Kt13(ru,:) = []; Kt13(:,ru) = [];
       % 4.
       dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
       Kt14 = Kts{matches(i)+4,2}; Kt14(ru,:) = []; Kt14(:,ru) = [];

       dksi = [0, 0, 0, 0, dksi11,dksi12,dksi13,dksi14];
       arclengths{i} = dksi;
       
       Ktprim0 = 1/epsil*(-1/2*Kt12 + 2*Kt11 - 3/2*Kt0);
       Ktprim11 = 1/epsil/2*(Kt12 - Kt0);
       Ktprim12 = 1/epsil/2*(Kt13 - Kt11);
       Ktprim13 = 1/epsil/2*(Kt14 - Kt12);
       Ktprim14 = 1/epsil*(3/2*Kt14 - 2*Kt13 + 1/2*Kt12);
       
       dKtprim0 = 1/epsil*(2*Kt0 - 5*Kt11 + 4*Kt12 - Kt13);

       [r0_,eigval0_] = solveCLEforMinEigNew(Kt0,Ktprim0,Kg,Kt0_0,typeofanal,matches(i));
       [r11_,eigval11_] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,typeofanal,matches(i)+1);
       [r12_,eigval12_] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,typeofanal,matches(i)+2);
       [r13_,eigval13_] = solveCLEforMinEigNew(Kt13,Ktprim13,Kg,Kt0_0,typeofanal,matches(i)+3);
       [r14_,eigval14_] = solveCLEforMinEigNew(Kt14,Ktprim14,Kg,Kt0_0,typeofanal,matches(i)+4);
           
       displacements_(:,5) = 0*Displ{matches(i)};
       displacements_(:,6) = Displ{matches(i)-0};
       displacements_(:,7) = Displ{matches(i)+1};
       displacements_(:,8) = Displ{matches(i)+2};
       displacements_(:,9) = Displ{matches(i)+3};
       displacements{i} = displacements_;
              
       if size(RR0,1)~=length(r0_)
           r0t = RR0; r11t = RR0; r12t = RR0; r13t = RR0; r14t = RR0;
           r0t(ra,:) = r0_;
           r11t(ra,:) = r11_;
           r12t(ra,:) = r12_;
           r13t(ra,:) = r13_;
           r14t(ra,:) = r14_;
       else
           r0t = r0_/norm(r0_);
           r11t = r11_/norm(r11_);
           r12t = r12_/norm(r12_);
           r13t = r13_/norm(r13_);
           r14t = r14_/norm(r14_);
       end
       
       
       for k = 1:size(r0_,2)
           if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
           if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
           if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end   
           if r13t(:,k)'*r14t(:,k)<0; r14t(:,k) = -r14t(:,k); end  
       end
       
       R = zeros(9,size(r0t,1),size(r0t,2));
       R(5,:,:) = r0t;
       R(6,:,:) = r11t;
       R(7,:,:) = r12t;
       R(8,:,:) = r13t;
       R(9,:,:) = r14t;
       
       EV = zeros(9,length(eigval0_));
       EV(5,:) = eigval0_;
       EV(6,:) = eigval11_;
       EV(7,:) = eigval12_;
       EV(8,:) = eigval13_;
       EV(9,:) = eigval14_;
       
    else
       % -4.
       Kt04 = Kts{matches(i)-4,2}; Kt04(ru,:) = []; Kt04(:,ru) = [];
       dksi04 = sqrt((Displ{matches(i)-4}-Displ{matches(i)-5})'*(Displ{matches(i)-4}-Displ{matches(i)-5}));
       % -3.
       Kt03 = Kts{matches(i)-3,2}; Kt03(ru,:) = []; Kt03(:,ru) = [];
       dksi03 = sqrt((Displ{matches(i)-3}-Displ{matches(i)-4})'*(Displ{matches(i)-3}-Displ{matches(i)-4}));
       % -2.
       Kt02 = Kts{matches(i)-2,2}; Kt02(ru,:) = []; Kt02(:,ru) = [];
       dksi02 = sqrt((Displ{matches(i)-2}-Displ{matches(i)-3})'*(Displ{matches(i)-2}-Displ{matches(i)-3}));
       % -1.
       Kt01 = Kts{matches(i)-1,2}; Kt01(ru,:) = []; Kt01(:,ru) = [];
       dksi01 = sqrt((Displ{matches(i)-1}-Displ{matches(i)-2})'*(Displ{matches(i)-1}-Displ{matches(i)-2}));
       % 0.
       Kt0 = Kts{matches(i),2}; Kt0(ru,:) = []; Kt0(:,ru) = [];
       % 1.
       dksi11 = sqrt((Displ{matches(i)}-Displ{matches(i)-1})'*(Displ{matches(i)}-Displ{matches(i)-1}));
       Kt11 = Kts{matches(i)+1,2}; Kt11(ru,:) = []; Kt11(:,ru) = [];
       % 2.
       dksi12 = sqrt((Displ{matches(i)+1}-Displ{matches(i)})'*(Displ{matches(i)+1}-Displ{matches(i)}));
       Kt12 = Kts{matches(i)+2,2}; Kt12(ru,:) = []; Kt12(:,ru) = [];
       % 3.
       dksi13 = sqrt((Displ{matches(i)+2}-Displ{matches(i)+1})'*(Displ{matches(i)+2}-Displ{matches(i)+1}));
       Kt13 = Kts{matches(i)+3,2}; Kt13(ru,:) = []; Kt13(:,ru) = [];
       % 4.
       dksi14 = sqrt((Displ{matches(i)+3}-Displ{matches(i)+2})'*(Displ{matches(i)+3}-Displ{matches(i)+2}));
       Kt14 = Kts{matches(i)+4,2}; Kt14(ru,:) = []; Kt14(:,ru) = [];
       
       dksi = [dksi04, dksi03, dksi02, dksi01, dksi11,dksi12,dksi13,dksi14];
       arclengths{i} = dksi;
       
       displacements_(:,1) = Displ{matches(i)-5};
       displacements_(:,2) = Displ{matches(i)-4};
       displacements_(:,3) = Displ{matches(i)-3};
       displacements_(:,4) = Displ{matches(i)-2};
       displacements_(:,5) = Displ{matches(i)-1};
       displacements_(:,6) = Displ{matches(i)-0};
       displacements_(:,7) = Displ{matches(i)+1};
       displacements_(:,8) = Displ{matches(i)+2};
       displacements_(:,9) = Displ{matches(i)+3};
       displacements{i} = displacements_;

       
       Ktprim03 = 1/(2*epsil)*(Kt02 - Kt04);
       Ktprim02 = 1/(2*epsil)*(Kt01 - Kt03);
       Ktprim01 = 1/(2*epsil)*(Kt0 - Kt02);
       Ktprim0 = 1/(2*epsil)*(Kt11 - Kt01);
       Ktprim11 = 1/(2*epsil)*(Kt12 - Kt0);
       Ktprim12 = 1/(2*epsil)*(Kt13 - Kt11);
       Ktprim13 = 1/(2*epsil)*(Kt14 - Kt12);
       
       dKtprim0 = 1/epsil*(Kt01 -2*Kt0 +Kt11);
       
       [r03_,eigval03_] = solveCLEforMinEigNew(Kt03,Ktprim03,Kg,Kt0_0,typeofanal,matches(i)-3);
       [r02_,eigval02_] = solveCLEforMinEigNew(Kt02,Ktprim02,Kg,Kt0_0,typeofanal,matches(i)-2);
       [r01_,eigval01_] = solveCLEforMinEigNew(Kt01,Ktprim01,Kg,Kt0_0,typeofanal,matches(i)-1);
       [r0_,eigval0_] = solveCLEforMinEigNew(Kt0,Ktprim0,Kg,Kt0_0,typeofanal,matches(i));
       [r11_,eigval11_] = solveCLEforMinEigNew(Kt11,Ktprim11,Kg,Kt0_0,typeofanal,matches(i)+1);
       [r12_,eigval12_] = solveCLEforMinEigNew(Kt12,Ktprim12,Kg,Kt0_0,typeofanal,matches(i)+2);
       [r13_,eigval13_] = solveCLEforMinEigNew(Kt13,Ktprim13,Kg,Kt0_0,typeofanal,matches(i)+3);
       
       
       if size(RR0,1)~=length(r0_)
           r03t = RR0; r02t = RR0; r01t = RR0; r0t = RR0; r11t = RR0; r12t = RR0; r13t = RR0;
           r03t(ra,:) = r03_;
           r02t(ra,:) = r02_;
           r01t(ra,:) = r01_;
           r0t(ra,:) = r0_;
           r11t(ra,:) = r11_;
           r12t(ra,:) = r12_;
           r13t(ra,:) = r13_;
       else
           r03t = r03_/norm(r03_);
           r02t = r02_/norm(r02_);
           r01t = r01_/norm(r01_);
           r0t = r0_/norm(r0_);
           r11t = r11_/norm(r11_);
           r12t = r12_/norm(r12_);
           r13t = r13_/norm(r13_);
       end
       
       
       for k = 1:size(r0_,2)
           if r03t(:,k)'*r02t(:,k)<0; r02t(:,k) = -r02t(:,k); end
           if r02t(:,k)'*r01t(:,k)<0; r01t(:,k) = -r01t(:,k); end
           if r01t(:,k)'*r0t(:,k)<0;  r0t(:,k)  = -r0t(:,k);  end
           if r0t(:,k)'*r11t(:,k)<0;  r11t(:,k) = -r11t(:,k); end
           if r11t(:,k)'*r12t(:,k)<0; r12t(:,k) = -r12t(:,k); end
           if r12t(:,k)'*r13t(:,k)<0; r13t(:,k) = -r13t(:,k); end   
       end
       
       R = zeros(9,size(r0t,1),size(r0t,2));
       R(2,:,:) = r03t;
       R(3,:,:) = r02t;
       R(4,:,:) = r01t;
       R(5,:,:) = r0t;
       R(6,:,:) = r11t;
       R(7,:,:) = r12t;
       R(8,:,:) = r13t;
       
       EV = zeros(9,length(eigval0_));
       EV(2,:) = eigval03_;
       EV(3,:) = eigval02_;
       EV(4,:) = eigval01_;
       EV(5,:) = eigval0_;
       EV(6,:) = eigval11_;
       EV(7,:) = eigval12_;
       EV(8,:) = eigval13_;
       
    end
    eigval{i} = EV;
    eigvec{i} = R;
    StiffMtxs{i,1} = Kt0;
    StiffMtxs{i,2} = Ktprim0;
    StiffMtxs{i,3} = dKtprim0;
    
%    [eq1,eq2] = solcheck(Kt0, typeofanal, Ktprim0, dKtprim0, lambda(matches(i)), epsil, EV, R, i, model);
%     for p = 1:10
%         r = R(5,:,p); r = r(:);
%         r(ru) = [];
%         L = EV(5,p);
%         disp(norm((Kt0 + L*Ktprim0)*r));
%     end
end

model.eigenvalues = eigval;
model.eigenvectors = eigvec;
model.energy = Energy;
model.arclengths = arclengths;
model.displacements = displacements;
model.lambda0 = lambda0';
model.stiffnessMatrices = StiffMtxs;

disp(['ready to save: ','AnalysisResults/',filename,'-',typeofanal,'.mat']);
save(['AnalysisResults/',filename,'-',typeofanal,'.mat'],'model');
disp('saved');
end

function [eq1,eq2] = solcheck(Kt, typeofanal, dKt, ddKt, lam, epsil, Ls, rs, i, model)
    eigpo = 1;
    
    L1 = Ls(5,eigpo) + lam;
    r01 = rs(4,:,eigpo); r01 = r01(:);  r01(model.BC(:,1)) = [];
    r0  = rs(5,:,eigpo); r0 = r0(:);    r0(model.BC(:,1)) = [];
    r11 = rs(6,:,eigpo); r11 = r11(:);  r11(model.BC(:,1)) = [];
    r12 = rs(7,:,eigpo); r12 = r12(:);  r12(model.BC(:,1)) = [];
    r13 = rs(8,:,eigpo); r13 = r13(:);  r13(model.BC(:,1)) = [];
    r1 = r0;
    if i == 1 % forward difference method
        dL1 = 1/epsil*(-3/2*Ls(5,eigpo) + 2*Ls(6,eigpo) - 1/2*Ls(7,eigpo));
        ddL1 = 1/epsil^2*(2*Ls(5,eigpo) - 5*Ls(6,eigpo) + 4*Ls(7,eigpo) - 1*Ls(8,eigpo));
        dr1 = 1/epsil*(-3/2*r0 + 2*r11 - 1/2*r12);
        ddr1 = 1/epsil^2*(2*r0 - 5*r11 + 4*r12 - 1*r13);
    else % central difference method
        dL1 = 1/epsil*(-1/2*Ls(4,eigpo) + 1/2*Ls(6,eigpo));
        ddL1 = 1/epsil^2*(1/2*Ls(4,eigpo) -2*Ls(5,eigpo)  + 1/2*Ls(6,eigpo));
        dr1 = 1/epsil*(-1/2*r01 + 1/2*r11);
        ddr1 = 1/epsil^2*(1/2*r01 -2*r0  + 1/2*r11);
    end


    switch typeofanal
        case 'I'
            B = -eye(size(Kt));
            dB = 0;
        case 'CLE'
            B = -dKt;
            dB = -ddKt;
        case 'I-CLE'
            B = eye(size(Kt)) - dKt;
            dB = -ddKt;
    end

    eq1   = norm((full(Kt) - (L1-lam)*B)*r1);
    eq2_1 = norm((dKt - (dL1-1)*B - (L1-lam)*dB)*r1);
    eq2_2 = norm((Kt - (L1-lam)*B)*dr1);
    eq2   = norm((dKt - (dL1-1)*B - (L1-lam)*dB)*r1 + (Kt - (L1-lam)*B)*dr1);

end
