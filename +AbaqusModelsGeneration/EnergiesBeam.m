function RES = EnergiesBeam(ELres_beam,beamNodes,beamElements)
if nargin<1
    name = 'testmodel';
    lambda = [0.1:0.1:1];
    geometry = [.3,.6,0.01,0.005,5];
    material = [210e9, 0.3];
    supports = {'C','C'};
    forces = [-100,20,2.0,50,5];
    [~, ~, beamNodes, beamElements] = Models.beamModel(name,lambda,geometry,material,supports,forces);
    [ELres_beam, ~] = AbaqusModelsGeneration.getHistoryOutputFromDatFile(['AbaqusRuns/',name,'_beam.dat']);
end

%% Logitudinal strain energy
strainIPs = ELres_beam('SE1');
strainIP1 = strainIPs(1:2:end,3:end);
strainIP2 = strainIPs(2:2:end,3:end);
stressIPs = ELres_beam('SF1');
stressIP1 = stressIPs(1:2:end,3:end);
stressIP2 = stressIPs(2:2:end,3:end);
EPSx_int1 = 0.5*strainIP1.*stressIP1;
EPSx_int2 = 0.5*strainIP2.*stressIP2;

zeromtx = zeros(size(EPSx_int1));

%% Shear xz strain energy
if ELres_beam.isKey('SE2')
    strainIPs = ELres_beam('SE2');
    strainIP1 = strainIPs(1:2:end,3:end);
    strainIP2 = strainIPs(2:2:end,3:end);
    stressIPs = ELres_beam('SF2');
    stressIP1 = stressIPs(1:2:end,3:end);
    stressIP2 = stressIPs(2:2:end,3:end);
    GAMMAxz_int1 = 0.5*strainIP1.*stressIP1;
    GAMMAxz_int2 = 0.5*strainIP2.*stressIP2;
else
    GAMMAxz_int1 = zeromtx;
    GAMMAxz_int2 = zeromtx;
end

%% Shear xy strain energy
if ELres_beam.isKey('SE3')
    strainIPs = ELres_beam('SE3');
    strainIP1 = strainIPs(1:2:end,3:end);
    strainIP2 = strainIPs(2:2:end,3:end);
    stressIPs = ELres_beam('SF3');
    stressIP1 = stressIPs(1:2:end,3:end);
    stressIP2 = stressIPs(2:2:end,3:end);
    GAMMAxy_int1 = 0.5*strainIP1.*stressIP1;
    GAMMAxy_int2 = 0.5*strainIP2.*stressIP2;
else
    GAMMAxy_int1 = zeromtx;
    GAMMAxy_int2 = zeromtx;
end

%% Bending y-y strain energy
strainIPs = ELres_beam('SK1');
strainIP1 = strainIPs(1:2:end,3:end);
strainIP2 = strainIPs(2:2:end,3:end);
stressIPs = ELres_beam('SM1');
stressIP1 = stressIPs(1:2:end,3:end);
stressIP2 = stressIPs(2:2:end,3:end);
XIy_int1 = 0.5*strainIP1.*stressIP1;
XIy_int2 = 0.5*strainIP2.*stressIP2;

%% Bending z-z strain energy
if ELres_beam.isKey('SK2')
    strainIPs = ELres_beam('SK2');
    strainIP1 = strainIPs(1:2:end,3:end);
    strainIP2 = strainIPs(2:2:end,3:end);
    stressIPs = ELres_beam('SM2');
    stressIP1 = stressIPs(1:2:end,3:end);
    stressIP2 = stressIPs(2:2:end,3:end);
    XIz_int1 = 0.5*strainIP1.*stressIP1;
    XIz_int2 = 0.5*strainIP2.*stressIP2;
else
    XIz_int1 = zeromtx;
    XIz_int2 = zeromtx;
end


%% Torsion (bending x-x) strain energy
if ELres_beam.isKey('SK3')
    strainIPs = ELres_beam('SK3');
    strainIP1 = strainIPs(1:2:end,3:end);
    strainIP2 = strainIPs(2:2:end,3:end);
    stressIPs = ELres_beam('SM3');
    stressIP1 = stressIPs(1:2:end,3:end);
    stressIP2 = stressIPs(2:2:end,3:end);
    XIx_int1 = 0.5*strainIP1.*stressIP1;
    XIx_int2 = 0.5*strainIP2.*stressIP2;
else
    XIx_int1 = zeromtx;
    XIx_int2 = zeromtx;
end

%% Bicurvature strain energy
if ELres_beam.isKey('BICURV')
    strainIPs = ELres_beam('BICURV');
    strainIP1 = strainIPs(1:2:end,3:end);
    strainIP2 = strainIPs(2:2:end,3:end);
    stressIPs = ELres_beam('BIMOM');
    stressIP1 = stressIPs(1:2:end,3:end);
    stressIP2 = stressIPs(2:2:end,3:end);
    biXI_int1 = 0.5*strainIP1.*stressIP1;
    biXI_int2 = 0.5*strainIP2.*stressIP2;
else
    biXI_int1 = zeromtx;
    biXI_int2 = zeromtx;
end
%%
epsx_energy = zeros(size(beamElements,1),size(strainIP1,2));
jxz_energy = epsx_energy; jxy_energy = epsx_energy;
xiy_energy = epsx_energy; xiz_energy = epsx_energy; xix_energy = epsx_energy;
bixi_energy = epsx_energy;
for i = 1:size(beamElements,1)
    elm = beamElements(i,:);
    n1coords = beamNodes(elm(2),2:end)';
    n2coords = beamNodes(elm(end),2:end)';
    
    epsx_int1 = EPSx_int1(i,:);
    epsx_int2 = EPSx_int2(i,:);
    Epsx = [epsx_int1', epsx_int2'];
    
    jxz_int1 = GAMMAxz_int1(i,:);
    jxz_int2 = GAMMAxz_int2(i,:);
    Jxz = [jxz_int1', jxz_int2'];
    
    jxy_int1 = GAMMAxy_int1(i,:);
    jxy_int2 = GAMMAxy_int2(i,:);
    Jxy = [jxy_int1', jxy_int2'];
    
    xiy_int1 = XIy_int1(i,:);
    xiy_int2 = XIy_int2(i,:);
    Xiy = [xiy_int1', xiy_int2'];
    
    xiz_int1 = XIz_int1(i,:);
    xiz_int2 = XIz_int2(i,:);
    Xiz = [xiz_int1', xiz_int2'];
    
    xix_int1 = XIx_int1(i,:);
    xix_int2 = XIx_int2(i,:);
    Xix = [xix_int1', xix_int2'];
    
    bixi_int1 = biXI_int1(i,:);
    bixi_int2 = biXI_int2(i,:);
    biXi = [bixi_int1', bixi_int2'];
    
    v1 = 0; v2 = sqrt(sum((n2coords-n1coords).^2));
    v = [v1;v2];
    for j = 1:size(Epsx,1)
        epsx_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Epsx(j,:));
        jxz_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Jxz(j,:));
        jxy_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Jxy(j,:));
        xiy_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Xiy(j,:));
        xiz_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Xiz(j,:));
        xix_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,Xix(j,:));
        bixi_energy(i,j) = AbaqusModelsGeneration.integrate1Dbeam(v,biXi(j,:));
    end

end
RES = containers.Map({'epsx','jxz','jxy','xiy','xiz','xix','bixi'},...
                     { epsx_energy , jxz_energy , jxy_energy ,...
                       xiy_energy , xiz_energy , xix_energy , bixi_energy });
end