function [optimised_parameters, fitness_values, optimised_parameters_history] = RLAPSO(NP, dim, range, x, y, z, max_iter, vel_clamp, w_min, w_max)
    % RLAPSO: 强化学习驱动的自适应粒子群优化算法（Reinforcement Learning-driven Adaptive Particle Swarm Optimization）
    % 输入参数：
    % NP              - 粒子数量
    % dim             - 参数维度
    % range           - 每个参数的最小值和最大值范围 (dim x 2)
    % x, y, z         - 输入数据（噪声数据）
    % max_iter        - 最大迭代次数
    % vel_clamp       - 速度限制因子
    % w_min, w_max    - 惯性权重的最小值和最大值

    epsilon = 0;  % 收敛阈值
    fitness_values = [];  % 用于记录每次迭代的全局最优适应度值
    optimised_parameters_history = zeros(max_iter, dim);  % 记录每轮迭代的参数

    % 获取参数范围
    min_range = range(:, 1);
    max_range = range(:, 2);

    % 粒子位置初始化
    sobol_seq = net(sobolset(dim), NP);
    x_pos = min_range' + (max_range - min_range)' .* sobol_seq;

    % 粒子速度初始化
    range_diff = (max_range - min_range);                            % dim x 1
    v_pos = vel_clamp * repmat(range_diff', NP, 1) .* rand(NP, dim); % NP x dim

    % 个人最优位置和适应度
    p_best_fitness = inf(1, NP);
    p_best = zeros(NP, dim);

    % 全局最优位置和适应度
    g_best_fitness = inf;
    g_best = zeros(1, dim);

    % 强化学习 Q-learning 初始化
    num_states = 10;  % 10量化搜索状态 5 8 12 15 20
    num_actions = 5;  % 动作数量
    Q_table = zeros(num_states, num_actions);
    alpha = 0.15;      % 学习率
    gamma = 0.95;      % 折扣因子
    %epsilon_rl = 0.2; % ε-greedy

    for itr = 1:max_iter
        % 计算每个粒子的适应度并更新个人最优和全局最优
        fitness = zeros(1, NP);
        for p = 1:NP
            param = x_pos(p, :);
            fitness(p) = objective_function(param, [x, y, z]);

            % 更新个人最优
            if fitness(p) < p_best_fitness(p)
                p_best_fitness(p) = fitness(p);
                p_best(p, :) = param;
            end

            % 更新全局最优
            if fitness(p) < g_best_fitness
                g_best_fitness = fitness(p);
                g_best = param;
            end
        end  

        % 自适应学习率
        beta = 1 / (0.05 * max_iter);  % 让学习率下降更慢
        alpha_t = alpha / (1 + beta * itr);

        % 自适应探索率
        epsilon_rl = 0.01 + (0.2 - 0.01) * exp(-0.05 * itr);

        % 动态调整惯性权重
        q1 = quantile(fitness, 0.25);   % 第1四分位数
        q3 = quantile(fitness, 0.75);   % 第3四分位数
        iqr = q3 - q1;                  % 四分位距
        outlier_threshold = q3 + 1.5 * iqr;  % 判定离群点的阈值
        outlier_ratio = sum(fitness > outlier_threshold) / NP;  % 离群点比例
        lambda = 3;  % 3控制增长曲率的参数 1 2 4 5
        w = w_min + (w_max - w_min) * (1 - exp(-lambda * outlier_ratio));

        % 强化学习：基于粒子的速度调整学习因子
        Psi_t = exp(-mean(abs(v_pos - circshift(v_pos, 1))));
        state_t = get_state_index(Psi_t, num_states);

        for p = 1:NP
            % 选择动作 (强化学习)
            action_t = select_action(Q_table, state_t, epsilon_rl);
            [delta_c1, delta_c2] = get_action_values(action_t);

            % 调整学习因子
            c1 = max(0.5, min(2.5, 1.5 + delta_c1(1)));  % 示例的初始值为 1.5
            c2 = max(0.5, min(2.5, 1.5 + delta_c2(1)));  % 示例的初始值为 1.5

            % 更新粒子速度
            v_pos(p, :) = w * v_pos(p, :) + ...
                              c1 * rand(1, dim) .* (p_best(p, :) - x_pos(p, :)) + ...
                              c2 * rand(1, dim) .* (g_best - x_pos(p, :));

            % 限制粒子速度范围
            v_pos(p, :) = max(min(v_pos(p, :), vel_clamp * (max_range - min_range)'), ...
                              -vel_clamp * (max_range - min_range)');

            % 更新粒子位置
            x_pos(p, :) = x_pos(p, :) + v_pos(p, :);
        end

        % 每 1 次迭代记录最优适应度值
        fitness_values = [fitness_values, g_best_fitness];

        if mod(itr, 10) == 0
            fprintf('迭代 %d / %d，当前最优适应度：%.4f\n', itr, max_iter, g_best_fitness);
        end

        optimised_parameters_history(itr, :) = g_best;  % 记录每轮全局最优参数

        % 收敛检测
        if g_best_fitness < epsilon
            fprintf('算法收敛于第 %d 次迭代\n', itr);
            break;
        end

        % 计算奖励值并更新 Q-table
        fit_improve = abs(g_best_fitness - min(fitness)) / (abs(min(fitness)) + eps);
        oscillation = mean(abs(v_pos(:) - circshift(v_pos(:), 1)));
        reward_t = compute_reward(fit_improve, oscillation);
        Q_table = update_q_table(Q_table, state_t, action_t, reward_t, alpha_t, gamma);
    end

    fprintf('达到最大迭代次数\n');
    optimised_parameters = g_best;
    fprintf('RLAPSO-拟合参数: \n');
    disp(optimised_parameters);
end

% 强化学习的辅助函数
function state_idx = get_state_index(Psi_t, num_states)
    Psi_t = mean(Psi_t(:));  % 确保 Psi_t 变成一个标量
    state_idx = min(num_states, max(1, round(Psi_t * num_states)));
end

function action_idx = select_action(Q_table, state_idx, epsilon)
    state_idx = state_idx(1);  % 取第一个数，确保是标量
    if rand < epsilon
        action_idx = randi(size(Q_table, 2));  % 随机选择
    else
        [~, action_idx] = max(Q_table(state_idx, :));  % 选择最优策略
    end
end

function [delta_c1, delta_c2] = get_action_values(action_idx)
    delta = 0.1;
    actions = [ delta, -delta;
               -delta,  delta;
                delta,  delta;
               -delta, -delta;
                0,      0 ];
    delta_c1 = actions(action_idx, 1);
    delta_c2 = actions(action_idx, 2);
end

function reward = compute_reward(fit_improve, oscillation)
    lambda1 = 1.5;  % 1.0适应度提升权重 0.5 1.0 1.0 1.5 1.5
    lambda2 = 1.5;    % 0.5震荡惩罚权重 0.5 1.0 1.5 1.0 1.5
    reward = lambda1 * fit_improve - lambda2 * oscillation;
end

function Q_table = update_q_table(Q_table, state_idx, action_idx, reward, alpha, gamma)
    max_future_q = max(Q_table(state_idx, :));
    Q_table(state_idx, action_idx) = ...
        (1 - alpha) * Q_table(state_idx, action_idx) + alpha * (reward + gamma * max_future_q);
end
