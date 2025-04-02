%Search for target keywords within files of a specific extension type under the target path.

fileList = dir(fullfile('D:\GraduateDissertation\Mang\MangAbaqus-master',  '**', '*.m'));  % 获取所有.m文件
for i = 1:length(fileList)
    filePath = fullfile(fileList(i).folder, fileList(i).name);  % 完整文件路径
    fileContent = fileread(filePath);  % 读取文件内容
    if contains(fileContent, 'u_max')   % 检查是否包含 "exec"
        fprintf('发现目标文件: %s\n', filePath);  % 高亮显示结果
    end
end