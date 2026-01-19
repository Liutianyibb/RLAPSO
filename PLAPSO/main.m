clc;
clear;
close all;

format longE

%% 生成模拟数据
% 设置模拟参数
% Hm = 0.5;                       % 地磁场强度 |M| = 0.5 高斯
% zaicHm = 49.705;                      % 数据集地磁场强度
Hm = 50000;                       % 地磁场强度 |M| = 50000 nT
D_angle = -5.8;                   % 磁偏角 
I_angle = 49.0;                   % 磁倾角
num_random_samples = 500;         % 生成1000个随机样本
noise_level_nT = 0;               % 选择噪声级别 (单位: nT, 可选0, 5, 10, 20, 50, 100)
outlier_percentage = 0.00;        % 选择离群数据比例

% 灵敏度误差
K1 = diag([1.0812, 0.9254, 1.0325]); 

% 非正交误差
% K2 = [   1,    0,   0.0028; 
%       0.0026,  1,   0.0035; 
%          0,    0,       1]; 

% alpha = 0.0231;
% beta  = 0.0512;
% gamma0 = 0.0378;
alpha = 0.075;
beta  = 0.057;
gamma0 = 0.123;

% K2 = [          cos(alpha),                      0,   sin(alpha);
%       sin(beta)*cos(gamma0),   cos(beta)*cos(gamma0),   sin(gamma0);
%                          0,                      0,           1];

K2 = [ cos(alpha)*cos(beta),    sin(alpha)*cos(beta),     sin(beta);
                          0,             cos(gamma0),   sin(gamma0);
                          0,                       0,            1];

K = K1 * K2;
fprintf('K:\n');
disp(K);

% 零偏误差
b = [101.26; 96.38; 110.54]; 

%调用函数生成加了干扰的数据
[x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan] = generate_magdata(num_random_samples, Hm, D_angle, I_angle, K1, K2, b, noise_level_nT, outlier_percentage);

%% 保存模拟数据
% % 创建表格，将数据和变量名对应
% data_table = table(x_noisy(:), y_noisy(:), z_noisy(:), x_yuan(:), y_yuan(:), z_yuan(:), ...
%                    'VariableNames', {'x_noisy', 'y_noisy', 'z_noisy', 'x_yuan', 'y_yuan', 'z_yuan'});
% 保存路径
save_path = 'D:\Download\Matlab2023b\xiangmu\DHEPSO\shuju\';
% 文件名
file_name1 = 'magdata(5nT).xlsx';
file_name2 = 'magdata(10nT).xlsx';
file_name3 = 'magdata(20nT).xlsx';
% 组合完整的文件路径
full_file_path = fullfile(save_path, file_name3);
% % 将表格数据写入 Excel 文件
% writetable(data_table, full_file_path, 'Sheet', 1, 'WriteVariableNames', true);


% % 绘制原始和带噪声椭球体数据
% figure;
% % 原始椭球体数据
% scatter3(x_yuan, y_yuan, z_yuan, 10, 'filled', 'r');    % 红色点
% hold on;
% % 带噪声和偏移的椭球体数据
% scatter3(x_noisy, y_noisy, z_noisy, 10, 'filled', 'b'); % 蓝色点
% % 设置图形属性
% title('原始与带噪声椭球体数据对比');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% legend('原始数据', '带噪声数据');
% axis equal;
% grid on;

%% 读取模拟数据
% % 读取 Excel 文件中的数据
% data_table = readtable(full_file_path, 'Sheet', 1);
% % 访问变量
% x_noisy = data_table.x_noisy;
% y_noisy = data_table.y_noisy;
% z_noisy = data_table.z_noisy;
% x_yuan = data_table.x_yuan;
% y_yuan = data_table.y_yuan;
% z_yuan = data_table.z_yuan;

%% 读取实际数据
% % 保存路径
% save_path0 = 'E:\论文\一种基于改进PSO算法和椭球拟合的自适应地磁信息校正方法\中国自动化大会\数据表\';
% % 文件名
% file_name0 = 'Data2.xlsx';
% % 组合完整的文件路径
% full_file_path0 = fullfile(save_path0, file_name0);
% % 读取 Excel 文件中的数据
% data_table = readtable(full_file_path0, 'Sheet', 1);
% % 访问变量
% x_noisy = data_table.x_noisy;
% y_noisy = data_table.y_noisy;
% z_noisy = data_table.z_noisy;

%% RLAPSO
 fitness_values_DHEPSO = RLAPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);
% fitness_values_RLAPSO = RLAPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, Hm);

%% PSO
%fitness_values_PSO = PSO_EllipsoidFit(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);
% fitness_values_PSO = PSO_EllipsoidFit(x_noisy, y_noisy, z_noisy, Hm);

%% MPSO
% fitness_values_MPSO = MPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);
% fitness_values_MPSO = MPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, Hm);

%% TSLPSO
% fitness_values_TSLPSO = TSLPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);

%% ADE
% fitness_values_ADE = ADE_EllipsoidFit(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);
% fitness_values_ADE = ADE_EllipsoidFit(x_noisy, y_noisy, z_noisy, Hm);

%% LSM
% [data_cali_LSM, data_raw_LSM] = LSM_EllipsoidFit(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm);
% 处理数据为列向量（确保按列保存）
fitness_col = fitness_values_DHEPSO(:); % 无论输入是行/列向量，均转换为列向量
save_path = 'DHEPSO_fitness.xlsx'; % 默认文件名
writematrix(fitness_col, save_path, 'Sheet', 1, 'Range', 'A1');
fprintf('适应度值已成功按列保存到：%s\n', fullfile(pwd, save_path));


