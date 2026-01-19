function [fitness_values_RLAPSO] = RLAPSO_Ellipsoid(x_noisy, y_noisy, z_noisy, x_yuan, y_yuan, z_yuan, Hm)

    %%%%%%%%%%%%%%%%%%%%%  RLAPSO-椭球拟合  %%%%%%%%%%%%%%%%%%%
    % 定义PSO参数
    NP = 100;             % 粒子数量
    dim = 9;              % 参数维度
    max_iter = 1000;      % 最大迭代次数
    v_clamp = 0.8;        % 速度限制因子
    w_min = 0.4;          % 0.4惯性权重的最小值 0.4 0.4 0.4 0.5 0.6 0.3
    w_max = 0.9;          % 0.9惯性权重的最大值 0.7 0.8 1.0 0.9 0.9 0.6

    % range = [-0.5, 0.5;      % 参数范围
    %          -0.5, 0.5; 
    %          -0.5, 0.5; 
    %           0.8, 1.2; 
    %           0.8, 1.2; 
    %           0.8, 1.2; 
    %          -500, 500; 
    %          -500, 500; 
    %          -500, 500];
    
    % range = [-1, 1;      % 参数范围
    %          -1, 1; 
    %          -1, 1; 
    %         0.5, 1.5; 
    %         0.5, 1.5; 
    %         0.5, 1.5; 
    %        -1000, 1000; 
    %        -1000, 1000; 
    %        -1000, 1000];

    % range = [-3, 3;      % 参数范围
    %          -3, 3; 
    %          -3, 3; 
    %         0.4, 1.6; 
    %         0.4, 1.6; 
    %         0.4, 1.6; 
    %       -1500, 1500; 
    %       -1500, 1500; 
    %       -1500, 1500];

    range = [-5, 5;      % 参数范围
             -5, 5; 
             -5, 5; 
            0.2, 1.8; 
            0.2, 1.8; 
            0.2, 1.8; 
          -2000, 2000; 
          -2000, 2000; 
          -2000, 2000];


    % 粒子群优化
    disp('开始粒子群优化...');
    tic; % 计时开始
    % 调用粒子群优化算法
    [optimised_parameters, fitness_values_RLAPSO, optimised_parameters_history] = RLAPSO(NP, dim, range, x_noisy, y_noisy, z_noisy, max_iter, v_clamp, w_min, w_max);
    % 显示耗时
    elapsed_time = toc; % 计时结束
    fprintf('优化耗时: %.2f 秒\n', elapsed_time);

    chazhi = optimised_parameters - [0.075, 0.057, 0.123, 1.0812, 0.9254, 1.0325, 101.26, 96.38, 110.54];
    fprintf('RLAPSO-拟合参数差值: \n');
    disp(chazhi);

    %% 数据校正
    % 提取拟合参数
    alpha = optimised_parameters(1);
    beta  = optimised_parameters(2);
    gamma0 = optimised_parameters(3);
    kx    = optimised_parameters(4);
    ky    = optimised_parameters(5);
    kz    = optimised_parameters(6);
    bx    = optimised_parameters(7);
    by    = optimised_parameters(8);
    bz    = optimised_parameters(9);

    % 构建非正交误差校准矩阵
    Kn = [ cos(alpha)*cos(beta),    sin(alpha)*cos(beta),     sin(beta);
                              0,             cos(gamma0),   sin(gamma0);
                              0,                       0,            1];

    % Kn = [   1,   0,   alpha;
    %       beta,   1,   gamma0;
    %          0,   0,      1];
    
    % 构建灵敏度误差矩阵
    Ksi = [kx,   0,    0;
            0,  ky,    0;
            0,   0,  kz];

    K = Ksi * Kn;
    
    % 构建零偏误差向量
    b = [bx; by; bz];

    % 数据校准
    % 初始化原始数据和校准后的数据
    data_raw_RLAPSO = [x_noisy, y_noisy, z_noisy]';
    data_yuan_RLAPSO = [x_yuan, y_yuan, z_yuan]';
    data_cali_RLAPSO = zeros(size(data_raw_RLAPSO));

    % 对每个数据点进行校正
    for i = 1:length(x_noisy)
        h = data_raw_RLAPSO(:,i);
        data_cali_RLAPSO(:,i) = inv(K) * (h - b);
    end

    % 输出校正参数
    fprintf('RLAPSO校正后K:\n');
    disp(K);
    fprintf('RLAPSO校正后H0:\n');
    disp(b);

    % 保存路径
    save_path = 'D:\Download\Matlab2023b\xiangmu\RLAPSO\shuju\';

    % % 原始数据保存为 Excel 文件
    % raw_data_file = fullfile(save_path, 'data_raw_RLAPSO(20nT).xlsx');    % 保存路径和文件名
    % writematrix(data_raw_RLAPSO', raw_data_file);  % 保存原始数据，转置为列格式
    % 
    % yuan_data_file = fullfile(save_path, 'data_yuan_RLAPSO(20nT).xlsx');    % 保存路径和文件名
    % writematrix(data_yuan_RLAPSO', yuan_data_file);  % 保存原始数据，转置为列格式
    % 
    % % 校正后数据保存为 Excel 文件
    % cali_data_file = fullfile(save_path, 'data_cali_RLAPSO(20nT).xlsx');    % 保存路径和文件名
    % writematrix(data_cali_RLAPSO', cali_data_file);  % 保存校正后数据，转置为列格式

    % % 提取各轴数据
    % x_cali = data_cali_RLAPSO(1,:)';
    % y_cali = data_cali_RLAPSO(2,:)';
    % z_cali = data_cali_RLAPSO(3,:)';
    % 
    % % 计算各轴的 MAE
    % mae_x = mean(abs(x_yuan - x_cali));
    % mae_y = mean(abs(y_yuan - y_cali));
    % mae_z = mean(abs(z_yuan - z_cali));
    % 
    % % 显示结果，以科学计数法保留两位小数
    % fprintf('x 轴的 MAE: %.4e\n', mae_x);
    % fprintf('y 轴的 MAE: %.4e\n', mae_y);
    % fprintf('z 轴的 MAE: %.4e\n', mae_z);

    %% 拟合结果绘制
    %%%%%%%%%%%%%%%  绘制DEL_PSO-椭球拟合结果图(十参数） %%%%%%%%%%%%
    % 绘制原始数据和校准后的数据的3D散点图
    figure
    scatter3(data_raw_RLAPSO(1,:), data_raw_RLAPSO(2,:), data_raw_RLAPSO(3,:), 10, 'blue', 'filled')        % 设置每个点的大小为10，原始数据，蓝色点
    hold on 
    scatter3(data_cali_RLAPSO(1,:), data_cali_RLAPSO(2,:), data_cali_RLAPSO(3,:), 10, 'red', 'filled')   % 设置每个点的大小为10，校正后的数据，红色点
    axis equal

    xlabel('X/uT')
    ylabel('Y/uT')
    zlabel('Z/uT')
    title('RLAPSO原始数据与校正后的数据');
    hold on 

    %% 创建一个新的图形窗口
    % figure;
    % 
    % subplot(3,1,1);
    % plot(data_cali_RLAPSO(1,:), 'r', 'LineWidth', 1.5);
    % grid on;
    % ylabel('Bx (\muT)');
    % title('校正后的 Bx 分量');
    % 
    % subplot(3,1,2);
    % plot(data_cali_RLAPSO(2,:), 'g', 'LineWidth', 1.5);
    % grid on;
    % ylabel('By (\muT)');
    % title('校正后的 By 分量');
    % 
    % subplot(3,1,3);
    % plot(data_cali_RLAPSO(3,:), 'b', 'LineWidth', 1.5);
    % grid on;
    % xlabel('样本点');
    % ylabel('Bz (\muT)');
    % title('校正后的 Bz 分量');
    % 
    % sgtitle('RLAPSO 校正后的磁场各分量变化曲线');


    %%%%%%%%% 绘制校正后的球 %%%%%%%%%
    % % 绘制校正后的球面
    % [u, v, w] = sphere(56);          % 生成一个由 (n+1)×(n+1) 个网格点构成的球面数据
    % K1 = Hm;                          % 校准后的磁场强度
    % shape1 = surf(K1*u, K1*v, K1*w);    % 根据拟合的半径K对球面进行缩放
    % set(shape1, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')  % 校正球体为红色，透明度为0.3
    % % shading flat                   % 设置表面渲染模式为平面着色
    % axis equal                       % 设置坐标轴单位长度相等

    %% 计算校正前后数据的模值
    raw_magnitudes = sqrt(sum(data_raw_RLAPSO.^2, 1));  % 原始数据的模值
    cali_magnitudes = sqrt(sum(data_cali_RLAPSO.^2, 1));  % 校正后数据的模值

    % 计算模值MAE
    res = 0;
    for i = 1:length(cali_magnitudes)

        magnitudes = cali_magnitudes(i);
        res = res + abs(magnitudes - 50000);

    end
    res_M = res / length(cali_magnitudes);
    disp(['RLAPSO模值平均绝对误差 (MAE): ', num2str(res_M)]);

    % % 保存模值数据
    % magnitudes_data_file = fullfile(save_path, 'magnitudes_data_RLAPSO(20nT).xlsx');    % 保存路径和文件名
    % 
    % % 创建包含标题和数据的表格
    % data_table = table((1:length(raw_magnitudes))', raw_magnitudes', cali_magnitudes', ...
    %                'VariableNames', {'索引', 'RLAPSO校正前', 'RLAPSO校正后'});

    % % 写入 Excel 文件
    % writetable(data_table, magnitudes_data_file);

    % 绘制模值对比图
    figure;
    plot(raw_magnitudes, 'r', 'DisplayName', '原始数据模值');
    hold on;
    plot(cali_magnitudes, 'b', 'DisplayName', '校正后数据模值');
    xlabel('数据点索引');
    ylabel('模值');
    legend show;
    title('RLAPSO校正前后数据的模值对比');
    grid on;

    %% 优化参数变化
    % figure;
    % for d = 1:dim
    %     subplot(ceil(dim/2), 2, d);  % 多子图排列
    %     plot(1:max_iter, optimised_parameters_history(:, d), 'LineWidth', 2);
    %     xlabel('迭代次数');
    %     ylabel(['参数 ', num2str(d)]);
    %     title(['参数 ', num2str(d), ' 变化趋势']);
    %     grid on;
    % end
    % sgtitle('优化参数随迭代变化曲线'); % 总标题



end
