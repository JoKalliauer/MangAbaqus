%#!
%university:TU Wien
%author: Modified for parallel computation of different h values
 close all 
 delete(findall(0,'type','figure','tag','TMWWaitbar'))
 set(0, 'DefaultFigureWindowState', 'normal');

  % there are following predefined test cases:
  %modelprops.testcase = 'TL_arch';
  modelprops.testcase = 'Kreis_arch3D'; % 使用Kreis_arch3D模型
  
  % possible element types (be aware of 2D and 3D):
  %3D
  % 定义需要分析的单元类型数组
  % %eltypes={'B32','B32H','B31','B31H','B33','B33H'}
  % % eltypes={'B32','B32H','B31','B33'} %B31H/B33H dofs aufpassen fuer rhoBungle
  %eltypes={'B32','B32H','B31','B33'}
  eltypes={'B32H'};
  %eltypes={'B32H'}
  
  % possible types of analysis
  % possible types of analysis
  %modelprops.typeofanalysis = 'I';modelprops.sigma=eps(1e-292); %identity matrix
  %modelprops.typeofanalysis = 'CLE';modelprops.sigma=pi() %
  %modelprops.typeofanalysis = 'KNL'; %[ (Kts+Ktu) - EW * Kt0 ] %konvergiert nicht
  modelprops.typeofanalysis = 'KNL2'; modelprops.sigma=0; %[ Kt - EW * Kt0 ]
  %modelprops.typeofanalysis = 'KNL3'; modelprops.sigma=1; %[ Kt0 + EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'KNL4'; modelprops.sigma=-1.1; %[ Kt0 - EW * (Kts+Ktu) ]
  %modelprops.typeofanalysis = 'Kg';
  %modelprops.typeofanalysisB = 'Kt0';
  %modelprops.typeofanalysisA = 'Ksigma';
  %modelprops.typeofanalysisA = 'KNoLinear';
  %modelprops.typeofanalysis=strcat(modelprops.typeofanalysisA,modelprops.typeofanalysisB);
  
  modelprops.numofelm = 10;
  
  epsil = 0.0005;%  0.01;
  sortType = 'none'; % eigenvectors sorting type: 'none', 'forwards', 'backwards'
  
  %plotfig= [14,28,33];
  %plotfig=[1:14,21,24,26,30,211];
  %plotfig=[2,7,14,21,26,211,30,34];
  %plotfig=[14,15,16,37,38,900,211];
  %plotfig=[7,14,15,30,211,43]; %#ok<NASGU>
  %plotfig=[15,943:945,948:949];main.savefigures=1
  %plotfig=[15,947,949,952,955:956,16,943,953,943,16];
  %plotfig=[15,45,35,19,52];%EW
  plotfig=[14,15,45,35];%EW
  forcedeig = []; %1; % forced eigenvector number
  
  modelprops.lambda = 0:epsil:1.5; % do not go over snap-through point
  modelprops.length=19.074;% [m] 
  
  modelprops.epsilon = epsil;
  modelprops.loadfactor = 1;
  
  modelprops.profil.tw = 8.6e-3;
  modelprops.forceAbaqus = 0; 
  modelprops.forcerun = 1; % false... do not force it; 0.5 force if it too less lambda, 1 ... always force it.
  modelprops.numofeigs = 18;
  modelprops.allowComplex = true;
  
  main.closall = false;
  main.savefigures = true;
  main.check = 0;
  main.colorshift = 0;
  modelprops.ask_delete = false;
  main.rstabil = 0.9999999960;
  modelprops.MeterValue = 1;
  main.whichEV = 'k0_11';
  main.Normierung = 'k0_11';
  main.rho = 'R1'; % KtR1 R1
  
  modelprops.followsigma = false;

  % 固定b值
  b_value = 10e-2 * modelprops.MeterValue;
  modelprops.profil.b = b_value;
  
  % 定义不同的h值数组
  %h_values = (10e-2:1e-2:30e-2) * modelprops.MeterValue;
  h_values = [27e-2:1e-2:29e-2, 30e-2:10e-2:100e-2] * modelprops.MeterValue;

  % 创建结果存储结构 - 为每种单元类型创建独立的存储空间
  results = cell(length(eltypes), length(h_values));
  models = cell(length(eltypes), length(h_values));
  
  % 检查是否支持并行计算
  try
      % 尝试创建一个并行池
      poolobj = gcp('nocreate');
      if isempty(poolobj)
          poolobj = parpool('local');
      end
      fprintf('检测到支持并行计算，将使用parfor循环进行计算\n');
      use_parfor = true;
  catch
      fprintf('未检测到并行计算支持或创建并行池失败，将使用普通for循环进行计算\n');
      use_parfor = false;
  end
  
  % 根据并行支持情况选择循环类型
  for et = 1:length(eltypes)
      % 设置当前单元类型
      modelprops.elementtype = char(eltypes(et));
      fprintf('开始计算单元类型 %s 的情况\n', eltypes{et});
      
      if use_parfor
          parfor i = 1:length(h_values)
              % 为每个并行任务创建独立的modelprops副本
              mp = modelprops;
              m = main;
              
              % 设置当前h值
              mp.profil.h = h_values(i);
              
              % 运行计算
              fprintf('计算单元类型 %s, h = %f的情况\n', eltypes{et}, h_values(i));
              [res, model] = Abaqus_single_run(mp, sortType, plotfig, forcedeig, m);
              
              % 存储结果
              results{et,i} = res;
              models{et,i} = model;
          end
      else
          for i = 1:length(h_values)
              % 为每个任务创建独立的modelprops副本
              mp = modelprops;
              m = main;
              
              % 设置当前h值
              mp.profil.h = h_values(i);
              
              % 运行计算
              fprintf('计算单元类型 %s, h = %f的情况\n', eltypes{et}, h_values(i));
              [res, model] = Abaqus_single_run(mp, sortType, plotfig, forcedeig, m);
              
              % 存储结果
              results{et,i} = res;
              models{et,i} = model;
          end
      end
  end

  
% 输出完成信息
fprintf('所有计算已完成\n');
fprintf('单元类型: ');
fprintf('%s ', eltypes{:});
fprintf('\nh值: ');
fprintf('%f ', h_values);
fprintf('\n结果已保存到Kreis_arch3D_results.mat\n');

% 创建保存目录
save_dir = 'D:/Data/MangAbaqus/Abaqus/Output/';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% 保存每个计算结果
for et = 1:length(eltypes)
    for i = 1:length(h_values)
        % 格式化文件名
        filename = sprintf('Kreis_arch3D_h%f_b%f_et%s_results.mat', h_values(i), b_value, eltypes{et});
        save_path = fullfile(save_dir, filename);
        
        % 保存当前参数组合的结果
        current_result = results{et,i};
        current_model = models{et,i};
        save(save_path, 'current_result', 'current_model', '-v7.3');
        fprintf('已保存结果到: %s\n', save_path);
    end
end