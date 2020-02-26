function RES = EnergiesShell(ELres_shell,shellNodes,shellElements)
if nargin<1
    name = 'testmodel';
    lambda = [0.1:0.1:1];
    geometry = [.3,.6,0.01,0.005,5];
    material = [210e9, 0.3];
    supports = {'C','C'};
    forces = [-100,20,2.0,50,5];
    [shellNodes, shellElements, ~, ~] = Models.beamModel(name,lambda,geometry,material,supports,forces);
    [ELres_shell, ~] = AbaqusModelsGeneration.getHistoryOutputFromDatFile(['AbaqusRuns/',name,'_shell.dat']);
end

%% Logitudinal x strain energy
strainIPs = ELres_shell('SE1');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF1');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
EPSx_int1 = 0.5*strainIP1.*stressIP1;
EPSx_int2 = 0.5*strainIP2.*stressIP2;
EPSx_int3 = 0.5*strainIP3.*stressIP3;
EPSx_int4 = 0.5*strainIP4.*stressIP4;

%% Logitudinal y strain energy
strainIPs = ELres_shell('SE2');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF2');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
EPSy_int1 = 0.5*strainIP1.*stressIP1;
EPSy_int2 = 0.5*strainIP2.*stressIP2;
EPSy_int3 = 0.5*strainIP3.*stressIP3;
EPSy_int4 = 0.5*strainIP4.*stressIP4;

%% Shear xy strain energy
strainIPs = ELres_shell('SE3');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF3');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
GAMMAxy_int1 = 0.5*strainIP1.*stressIP1;
GAMMAxy_int2 = 0.5*strainIP2.*stressIP2;
GAMMAxy_int3 = 0.5*strainIP3.*stressIP3;
GAMMAxy_int4 = 0.5*strainIP4.*stressIP4;

%% Shear xz strain energy
strainIPs = ELres_shell('SE4');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF4');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
GAMMAxz_int1 = 0.5*strainIP1.*stressIP1;
GAMMAxz_int2 = 0.5*strainIP2.*stressIP2;
GAMMAxz_int3 = 0.5*strainIP3.*stressIP3;
GAMMAxz_int4 = 0.5*strainIP4.*stressIP4;

%% Shear yz strain energy
strainIPs = ELres_shell('SE5');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF5');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
GAMMAyz_int1 = 0.5*strainIP1.*stressIP1;
GAMMAyz_int2 = 0.5*strainIP2.*stressIP2;
GAMMAyz_int3 = 0.5*strainIP3.*stressIP3;
GAMMAyz_int4 = 0.5*strainIP4.*stressIP4;

%% Through deepth zz strain energy
strainIPs = ELres_shell('SE6');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SF6');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
EPSz_int1 = 0.5*strainIP1.*stressIP1;
EPSz_int2 = 0.5*strainIP2.*stressIP2;
EPSz_int3 = 0.5*strainIP3.*stressIP3;
EPSz_int4 = 0.5*strainIP4.*stressIP4;

%% Bending x-x strain energy
strainIPs = ELres_shell('SK1');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SM1');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
XIx_int1 = 0.5*strainIP1.*stressIP1;
XIx_int2 = 0.5*strainIP2.*stressIP2;
XIx_int3 = 0.5*strainIP3.*stressIP3;
XIx_int4 = 0.5*strainIP4.*stressIP4;

%% Bending y-y strain energy
strainIPs = ELres_shell('SK2');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SM2');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
XIy_int1 = 0.5*strainIP1.*stressIP1;
XIy_int2 = 0.5*strainIP2.*stressIP2;
XIy_int3 = 0.5*strainIP3.*stressIP3;
XIy_int4 = 0.5*strainIP4.*stressIP4;

%% Bending xy strain energy
strainIPs = ELres_shell('SK3');
strainIP1 = strainIPs(1:4:end,3:end);
strainIP2 = strainIPs(2:4:end,3:end);
strainIP3 = strainIPs(3:4:end,3:end);
strainIP4 = strainIPs(4:4:end,3:end);
stressIPs = ELres_shell('SM3');
stressIP1 = stressIPs(1:4:end,3:end);
stressIP2 = stressIPs(2:4:end,3:end);
stressIP3 = stressIPs(3:4:end,3:end);
stressIP4 = stressIPs(4:4:end,3:end);
XIxy_int1 = 0.5*strainIP1.*stressIP1;
XIxy_int2 = 0.5*strainIP2.*stressIP2;
XIxy_int3 = 0.5*strainIP3.*stressIP3;
XIxy_int4 = 0.5*strainIP4.*stressIP4;

%%
epsx_energy = zeros(size(shellElements,1),size(strainIP1,2));
epsy_energy = epsx_energy; epsz_energy = epsx_energy;
jxy_energy = epsx_energy; jxz_energy = epsx_energy; jyz_energy = epsx_energy;
xix_energy = epsx_energy; xiy_energy = epsx_energy; xixy_energy = epsx_energy; 

for i = 1:size(shellElements,1)
    elm = shellElements(i,:);
    n1coords = shellNodes(elm(2),2:end)';
    n2coords = shellNodes(elm(3),2:end)';
    n3coords = shellNodes(elm(4),2:end)';
    n4coords = shellNodes(elm(5),2:end)';
    
    xVect = [1, 0, 0]; yVect = [0, 1, 0]; zVect = [0, 0, 1]; 
    normalVect = cross((n3coords - n1coords),(n4coords - n2coords));
    zVect_ = normalVect/norm(normalVect);
    xVect_ = [1, 0, 0]';
    yVect_ = cross(normalVect,xVect); yVect_ = yVect_'/norm(yVect_);
    
    T = [xVect*xVect_, yVect*xVect_, zVect*xVect_;
         xVect*yVect_, yVect*yVect_, zVect*yVect_;
         xVect*zVect_, yVect*zVect_, zVect*zVect_];
%     if zVect_(3)~=1
%         disp('tu')
%     end
    n1coords = T*n1coords; n1coords = n1coords(1:2);
    n2coords = T*n2coords; n2coords = n2coords(1:2);
    n3coords = T*n3coords; n3coords = n3coords(1:2);
    n4coords = T*n4coords; n4coords = n4coords(1:2);
    
    epsx_int1 = EPSx_int1(i,:);
    epsx_int2 = EPSx_int2(i,:);
    epsx_int3 = EPSx_int3(i,:);
    epsx_int4 = EPSx_int4(i,:);
    Epsx = [epsx_int1', epsx_int2', epsx_int3', epsx_int4'];
    
    epsy_int1 = EPSy_int1(i,:);
    epsy_int2 = EPSy_int2(i,:);
    epsy_int3 = EPSy_int3(i,:);
    epsy_int4 = EPSy_int4(i,:);
    Epsy = [epsy_int1', epsy_int2', epsy_int3', epsy_int4'];
    
    epsz_int1 = EPSz_int1(i,:);
    epsz_int2 = EPSz_int2(i,:);
    epsz_int3 = EPSz_int3(i,:);
    epsz_int4 = EPSz_int4(i,:);
    Epsz = [epsz_int1', epsz_int2', epsz_int3', epsz_int4'];
    
    jxz_int1 = GAMMAxz_int1(i,:);
    jxz_int2 = GAMMAxz_int2(i,:);
    jxz_int3 = GAMMAxz_int3(i,:);
    jxz_int4 = GAMMAxz_int4(i,:);
    Jxz = [jxz_int1', jxz_int2', jxz_int3', jxz_int4'];
    
    jxy_int1 = GAMMAxy_int1(i,:);
    jxy_int2 = GAMMAxy_int2(i,:);
    jxy_int3 = GAMMAxy_int3(i,:);
    jxy_int4 = GAMMAxy_int4(i,:);
    Jxy = [jxy_int1', jxy_int2', jxy_int3', jxy_int4'];
    
    jyz_int1 = GAMMAyz_int1(i,:);
    jyz_int2 = GAMMAyz_int2(i,:);
    jyz_int3 = GAMMAyz_int3(i,:);
    jyz_int4 = GAMMAyz_int4(i,:);
    Jyz = [jyz_int1', jyz_int2', jyz_int3', jyz_int4'];
    
    xix_int1 = XIx_int1(i,:);
    xix_int2 = XIx_int2(i,:);
    xix_int3 = XIx_int3(i,:);
    xix_int4 = XIx_int4(i,:);
    Xix = [xix_int1', xix_int2', xix_int3', xix_int4'];
    
    xiy_int1 = XIy_int1(i,:);
    xiy_int2 = XIy_int2(i,:);
    xiy_int3 = XIy_int3(i,:);
    xiy_int4 = XIy_int4(i,:);
    Xiy = [xiy_int1', xiy_int2', xiy_int3', xiy_int4'];
    
    xixy_int1 = XIxy_int1(i,:);
    xixy_int2 = XIxy_int2(i,:);
    xixy_int3 = XIxy_int3(i,:);
    xixy_int4 = XIxy_int4(i,:);
    Xixy = [xixy_int1', xixy_int2', xixy_int3', xixy_int4'];

    v = [n1coords'; n2coords'; n3coords'; n4coords'];
    for j = 1:size(Epsx,1)
        epsx_energy(i,j) = integrate2Dquad(v,Epsx(j,:));
        epsy_energy(i,j) = integrate2Dquad(v,Epsy(j,:));
        epsz_energy(i,j) = integrate2Dquad(v,Epsz(j,:));
        
        jxy_energy(i,j) = integrate2Dquad(v,Jxy(j,:));
        jxz_energy(i,j) = integrate2Dquad(v,Jxz(j,:));
        jyz_energy(i,j) = integrate2Dquad(v,Jyz(j,:));
        
        xix_energy(i,j) = integrate2Dquad(v,Xix(j,:));
        xiy_energy(i,j) = integrate2Dquad(v,Xiy(j,:));
        xixy_energy(i,j) = integrate2Dquad(v,Xixy(j,:));
    end

end
RES = containers.Map({'epsx','epsy','epsz','jxy','jxz','jyz','xix','xiy','xixy'},...
                     { epsx_energy , epsy_energy , epsz_energy ,...
                       jxy_energy, jxz_energy , jyz_energy ,...
                       xix_energy, xiy_energy, xixy_energy});
end