%% main
% 初始化
map = 0; 
derta = 15; % map size = 15*15
start_x = 100;
start_y = 100;
generate_map(map, start_x, start_y, derta, derta); 
load(sprintf('output/map_%d/map%d.mat', map, map));

% 设置 sample_rate 从 0.01 到 0.15，生成 15 个等间距的数值
sample_rate_values = linspace(0.01, 0.15, 15);

% 生成 psi_random 和 psi_mmi
psi_file = sprintf('output/map_%d/psi%d.mat', map, map);

if exist(psi_file, 'file')  % 如果文件存在，直接加载
    fprintf('psi already exists.\n');
    load(psi_file, 'psi', 'psi_opt');
else  % 如果文件不存在，生成并保存
    [psi, psi_opt] = generate_psi(sample_rate_values, M, N, phi);  % 注意：不再传递 psi_file 参数
    save(psi_file, 'psi', 'psi_opt');  % 保存生成的 psi 和 psi_opt
end

% 选择是否进行恢复算法
result_file = sprintf('output/map_%d/result%d.mat', map, map);
if ~exist(result_file, 'file')
    % 设置参数
    sigma = 0.1; % Noise power
    noise = randn(M, 1) * sigma; % Noise component
    
    % 设置迭代次数
    iteration_time_omp = 100;
    iteration_time_sbl = 50;
    
    % 创建输出结构体
    result_random_omp = struct('omega_est', [], 'mse', [],'converge_point', [], 'rn_norm', []);
    result_random_sbl = struct('omega_est', [], 'mse', [], 'beta_cur', []);
    result_mmi_sbl = struct('omega_est', [], 'mse', [], 'beta_cur', []);
    
    % 遍历不同的 sample_rate 值并计算相应的 MSE
    for i = 1:length(sample_rate_values)
        
        % Random_OMP
        [omega_est, mse, rn_norm, converge_point] = recover_omp(phi, psi{i}, omega_real, noise, iteration_time_omp);
        result_random_omp(i) = struct('omega_est', omega_est, 'mse', mse, 'rn_norm', rn_norm, 'converge_point', converge_point); 
    
        % Random_SBL
        [omega_est, mse, beta_cur] = recover_sbl(phi, psi{i}, omega_real, noise, iteration_time_sbl, sigma);
        result_random_sbl(i) = struct('omega_est', omega_est, 'mse', mse, 'beta_cur', beta_cur);
    
        % MMI_SBL
        [omega_est, mse, beta_cur] = recover_sbl(phi, psi_opt{i}, omega_real, noise, iteration_time_sbl, sigma); 
        result_mmi_sbl(i) = struct('omega_est', omega_est, 'mse', mse, 'beta_cur', beta_cur);
    
    end

    save(sprintf('output/map_%d/result%d.mat', map, map), ...
        'result_random_omp','result_random_sbl','result_mmi_sbl','sample_rate_values','sigma','noise');

else
    load(result_file);
end


% Plot part
command = input('Plot or not?\n', 's');
if strcmpi(command, 'yes')

    plot_basic(map); % 地图基本信息：RSS_entire / RSS_roi
    plot_mse(map); % 恢复信号 MSE vs Sampling_rate
    
    for sample_rate_set = 1:15 % 选择采样率 0.01 ~ 0.15
        rate = sample_rate_set*0.01;
        img_dir = sprintf('output/map_%d/sample_rate=%.2f', map, rate);
        if ~exist(img_dir,'dir')
            mkdir(img_dir);
        end
        
        % 恢复信号图
        plot_recover_signal(map, sample_rate_set);
        saveas(gcf, fullfile(img_dir, sprintf('RandomOMP_rate=%.2f.png', rate)));
    
        % OMP收敛性图
        plot_converge_omp(map, sample_rate_set);
        saveas(gcf, fullfile(img_dir, sprintf('RandomSBL_rate=%.2f.png', rate)));
    
        % SBL收敛性图
        plot_converge_sbl(map, sample_rate_set);
        saveas(gcf, fullfile(img_dir, sprintf('Recovered_Map-rate=%.2f.png', rate)));
    
    end

else
    return;
end



