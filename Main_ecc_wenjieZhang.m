%#!
%university:TU Wien
%author of this script: Johannes Kalliauer(2020-2025)
%author of subprograms: Johannes Kalliauer(2020-2025), Michal Malendowski (2019-2020)
%created: ~2020

%% Last changes
%2025-01-xx modelprops.loadfactor = 0;

%% Code

 %#ok<*NOPTS>
 %clear
 close all
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
  
 
  %[~,modelprops.ecc]=eccfromU(0.5)
  %modelprops.ecc=0.0403316;
  modelprops.testcase = 'ec2CompressionBeam'; 

  
  modelprops.length = 5;
  
  eltypes={'B32OSH'};
 

  modelprops.typeofanalysis = 'KNL2'; %modelprops.sigma=0;
 
  sortType = 'none';  %sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards','forwardJK'
  %plotfig=[2,14,15,19,35,45,976,977];
  plotfig=[1:15,19,52];
  forcedeig = []; %1; % forced eigenvector number 'none' sorting

  
  modelprops.loadfactor = 1;
  
  modelprops.profil.tw= 8.6e-3;
  modelprops.forceAbaqus=0; %-1..returns error if not exist, 0..use old if exist, 1.. force new calc
  modelprops.forcerun=1; %0..use existing one, 0.5.. force run if last lambda smaller than requested, always fore a new calc.
  modelprops.allowComplex=1;%0..no complex, 1 also complex, 2 only complex
  main.closall=true;
  main.savefigures=1;
  main.check=0;
  main.colorshift=0;
  modelprops.ask_delete=true;
  main.rsame=NaN;%0.8;
  main.rstabil=NaN;%0.99999;
  main.whichEV='k0_11'; % main.whichEV='bungle'; 'Disp'; 'Rot'; 'wrap'; 'Hyb'; 'skip' ; 'k11' 'k0_11' '2023-12' '2023_12Hyb' '2023_12noHyb' '2023_12half'
  main.Normierung='k0_11'; % 'R1';  'rNCT_K0_r' ; 'rCT_K0_r' 'sqrtK_r' 'skip' 'k0_11'
  main.rho='R1'; % KtR1 R1 'A0R1' 
  
  %modelprops.MeterValue=1; %1000mm=1m=0.001km ;
  main.xBezug='n'; %n..normalisiert; d..differenz zut Refwert; 1...Abaqus-Lambda; s...Stepnumber; i..individual
  main.flipAxis=false;
  
  modelprops.sigma=0;
  modelprops.followsigma=1;
  modelprops.sortJKeigval=-1; %1..closest to zero, -1 ..most negative one
  
numofelms={2};


%Exz={modelprops.ecc};
modelprops.numofeigs=7;%min 7 EV

%epsils={0.005}% nicht kleiner als 0.0005

%U_values = 0:0.05:0.9999;%到0.95
U_values =0.96:0.01:0.99;
modelprops.ecc = zeros(1, numel(U_values));
epsils_list={0.001}
for k = 1:numel(U_values)
    [~,tmp] = eccfromU(U_values(k));  
    modelprops.ecc(k) = tmp; 
end

% 初始化并行计算环境
if ~exist('gcp','file')
    % 如果Parallel Computing Toolbox不可用，使用普通for循环
    useParallel = false;
else
    % 尝试启动并行池
    try
        parpool('local');
        useParallel = true;
    catch
        useParallel = false;
        warning('无法启动并行池，将使用普通for循环');
    end
end

if useParallel
    % 使用parfor并行计算
    parfor n=1:numel(epsils_list)
        % 创建局部变量以避免并行计算中的冲突
        local_modelprops = modelprops;
        local_main = main;
        local_epsils = epsils_list(n);
        
        for num=1:numel(U_values)
            local_Exz = {local_modelprops.ecc(num)};
            for l=1:numel(local_epsils)
                local_modelprops.epsilon = cell2mat(local_epsils(l));
                local_modelprops.lambda = 0:local_modelprops.epsilon:max([2.5,20*local_modelprops.epsilon]);
                for k=1:numel(local_Exz)
                    for j=1:numel(numofelms)
                        local_modelprops.numofelm = cell2mat(numofelms(j));
                        for i=1:numel(eltypes)
                            local_elementtype = char(eltypes(i));
                            local_main.colorshift=(i-1)+(j-1)*numel(eltypes)+(l-1)*numel(eltypes)*numel(numofelms);
                            [~,~] = Abaqus_single_run(local_modelprops,sortType,plotfig,forcedeig,local_main,local_modelprops.numofelm,cell2mat(local_Exz(k)),local_elementtype);
                        end
                    end
                end
            end
        end
    end
else
    % 使用普通for循环
    for n=1:numel(epsils_list)
        epsils=epsils_list(n);
        for num=1:numel(U_values)
            Exz={modelprops.ecc(num)};
            for l=1:numel(epsils)
                modelprops.epsilon = cell2mat(epsils(l));
                modelprops.lambda = 0:modelprops.epsilon:max([2.5,20*modelprops.epsilon]);
                for k=1:numel(Exz)
                    for j=1:numel(numofelms)
                        modelprops.numofelm = cell2mat(numofelms(j));
                        for i=1:numel(eltypes)
                            elementtype = char(eltypes(i));
                            main.colorshift=(i-1)+(j-1)*numel(eltypes)+(l-1)*numel(eltypes)*numel(numofelms);
                            [res,model] = Abaqus_single_run(modelprops,sortType,plotfig,forcedeig,main,modelprops.numofelm,cell2mat(Exz(k)),elementtype);
                        end
                    end
                end
            end
        end
    end
end

allFigures = findall(0,'Type','figure'); % find all figures
set(allFigures,'WindowState','normal'); % set the WindowState of all figures to normal

