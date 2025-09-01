%% main
clear; close all; clc;

%% Set params
mapname         = "tinymap"    ; % {tinymap, largemap, hugemap}
K               = 4             ; % {4, 8, 12, 16}
sample_rate     = 0.15          ; % (sa) 0.01 ~ 0.15 / 0.05, 0.10, 0.15
select_rate     = 0.5           ; % (se) Select partial Rx
method_phi      = "idw"         ; % {idw, halrtc, kriging}
method_psi      = "random"      ; % {random, mmi}
method_recov    = "cmsbl"         ; % {omp, sbl, csbl, msbl, cmsbl}
sigma2          = 0.05          ; % Noise power: σ^2

% Output directory
if ~exist("result", "dir")
    mkdir("result");
end

dir = sprintf("result/%s", mapname);
if ~exist(dir,"dir")
    mkdir(dir);
end

% Code start
fprintf('---- Code started ----\n');
method = sprintf('%s_%s_%s', method_phi, method_psi, method_recov); % entire method name
exp = sprintf('%s_%s_K=%d_sa=%.2f_se=%.2f_si=%.2f', mapname, method, K, sample_rate, select_rate, sigma2);
fprintf('Map: %s\n', mapname);
fprintf('Method: %s\n', method);

% Create output directory
dir_out = sprintf("%s/%s", dir, exp);
if ~exist(dir_out,"dir")
    mkdir(dir_out);
end

% Load map file
file_map = sprintf('%s/%s.mat', dir, mapname);
if ~exist(file_map, "file")
    fprintf("Generating map: %s\n", mapname);
    map = generate_map2D(mapname);
    save(file_map, 'map'); 
else
    fprintf("Read map file: %s\n", file_map);
    load(file_map);
end

% Load map file with specific selected rate
file_map_se = sprintf('%s/%s_K=%d_se=%.2f.mat', dir, mapname, K, select_rate);
if ~exist(file_map_se, "file")
    fprintf("Selected rate = %.2f in map: %s\n",select_rate, mapname);

    % Select K Tx
    map.sparsity   = K;
    map.select_rate = select_rate;
    
    % Select partial Rx
    map.validPos = setdiff(1:map.size, map.build)';
    map.validPosNum = numel(map.validPos);
    map.partPos = map.validPos(randperm(map.validPosNum, floor(select_rate * map.validPosNum)));
    map.interPos = setdiff(map.validPos, map.partPos);
    
    % Construct ω_real
    map.selectedTxPos = [map.Tx(1:K).pos];
    map.omega_real = sparse(map.selectedTxPos, 1, [map.Tx(1:K).Pt_mW], map.size, 1);

    save(file_map_se, 'map');
else
    fprintf("Read map file: %s\n", file_map_se);
    load(file_map_se);
end

%% Generate sparse dictionary
file_phi = sprintf('%s/%s_phi_%s_K=%d_se=%.2f.mat', dir, ...
    mapname, method_phi, K, select_rate);
if ~exist(file_phi, "file")
    fprintf("Generating φ by %s\n", method_phi);
    [phi, phi_rt, phi0] = generate_phi(method_phi, map);
    save(file_phi,'method_phi', 'phi', 'phi_rt', 'phi0');
else 
    fprintf("Read φ file: %s\n", file_phi);
    load(file_phi);
end

%% Generate measurement matrix
file_psi = sprintf('%s/%s_psi_%s_K=%d_sa=%.2f_se=%.2f.mat', dir_out, ...
    mapname, method_psi, K, sample_rate, select_rate);
if ~exist(file_psi, "file")
    fprintf("Generating φ by %s\n", method_psi);
    [psi] = generate_psi(method_psi, map, sample_rate, phi);
    save(file_psi,'method_psi', 'psi');
else 
    fprintf("Read ψ file: %s\n", file_psi);
    load(file_psi);
end

% Transmit process
omega_real = map.omega_real;
Phi = psi * phi; % Sensing matrix
y = Phi * omega_real; % Observation vector

%% Recover signal
[omega_est] = recover_signal(method_recov, y, Phi, sigma2, K);

rss_full_vec = phi0 * map.omega_real;
rss_real_vec = phi_rt * map.omega_real;
rss_est_vec = phi * omega_est;

rss_full_dbm = ln_to_db(rss_full_vec);
rss_real_dbm = ln_to_db(rss_real_vec);
rss_est_dbm = ln_to_db(rss_est_vec);

rss_full2D = reshape(rss_full_dbm, map.Nx, map.Ny);
rss_real2D = reshape(rss_real_dbm, map.Nx, map.Ny);
rss_est2D = reshape(rss_est_dbm, map.Nx, map.Ny);

%% Evaluation
valid.mask = isfinite(rss_est_dbm) & isfinite(rss_full_dbm); % ignore Inf value
valid.num = numel(valid.mask);

% MSE - sparse signal recovery
mse = 10 * log10(norm(omega_real - omega_est) / norm(omega_real));
% RMSE - RSS recovery
rmse = norm(rss_est_dbm(valid.mask) - rss_full_dbm(valid.mask)) / sqrt(valid.num);
% MAE - dictionary
mae = norm(rss_est_dbm(valid.mask) - rss_full_dbm(valid.mask), 1) / valid.num;

fprintf('\n---- Evaluation: %s ----\n', exp);
fprintf('MSE = %.4f dB\n', mse);
fprintf('RMSE = %.4f dB\n', rmse);
fprintf('MAE = %.4f dB\n', mae);
clear("valid");

%% Plot and save figure
figName = {
    sprintf('%s-Ground_Truth.png', map.name), 
    sprintf('%s_K=%d_se=%.2f.png', map.name, K, select_rate), 
    sprintf('%s.png', exp)
};
plotRSSfig(map, rss_full2D, rss_real2D, rss_est2D, figName, dir);
clear("figName");

% Save result
result_name = sprintf('%s/%s.mat', dir_out, exp);
save(result_name);

fprintf('\n---- Code finished ----\n');

%% Additional function
% Convert dB to linear
function ln = db_to_ln(db)
    ln = 10.^(db./10);
end

% Convert linear to dB
function db = ln_to_db(ln)
    db = 10.*log10(ln);
end

% Convert index to coordinate - 2D
function [x, y] = itc(idx, Nx)
    x = mod(idx-1, Nx) + 1;
    y = floor((idx-1) / Nx) + 1;
end

% Plot and save RSS figure
function plotRSSfig(map, rss_full2D, rss_real2D, rss_est2D, figName, dir)
    % 初始化设置
    [buildX, buildY] = itc(map.build, map.Nx);
    n = 512;
    h = linspace(0.7, 0, n)';
    s = linspace(0.8, 1, n)';
    v = linspace(0.9, 1, n)';
    cmap = hsv2rgb([h, s, v]);
    barMin = -75;
    barMax = 10;
    
    % 创建主图窗
    fig = figure('Visible', 'off', 'Position', [250, 100, 1200, 700]);

    % 检查文件是否存在
    findFig = [
        exist(fullfile(dir, figName{1}), 'file'), 
        exist(fullfile(dir, figName{2}), 'file'),
        exist(fullfile(dir, figName{3}), 'file')
    ];
    
    % 遍历三个子图
    titles = {'Ground Truth', 'RSS real', 'RSS est'};
    data = {rss_full2D, rss_real2D, rss_est2D};
    
    for k = 1:3
        % 创建子图并绘制数据
        ax = subplot(1,3,k);
        imagesc(data{k});
        set(ax, 'YDir','normal');
        axis([1 map.Nx 1 map.Ny]);
        axis equal tight;
        hold(ax, 'on');
        colormap(ax, cmap);
        set(ax, 'CLim', [barMin barMax]); % 统一颜色范围
        
        % 添加Colorbar并获取句柄
        cb = colorbar(ax, 'Location', 'eastoutside'); % 强制右侧显示
        ylabel(cb, 'dBm');
        title(titles{k});
        ylabel('X axis');
        
        % 计算散点尺寸
        dpi = get(0, 'ScreenPixelsPerInch');
        figPos = get(gcf, 'Position');
        axPos = ax.Position;
        axWidthInches = figPos(3) * axPos(3) / dpi;
        axHeightInches = figPos(4) * axPos(4) / dpi;
        cellWidthInches = axWidthInches / map.Nx;
        cellHeightInches = axHeightInches / map.Ny;
        markerDiameterInches = max(cellWidthInches, cellHeightInches);
        markerAreaPointsSquared = (markerDiameterInches * 72)^2 * 0.98;
        
        % 绘制散点
        scatter(ax, buildY, buildX, markerAreaPointsSquared, 's', ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
        
        % 同步Colorbar与图像高度
        axPos = ax.Position; % 获取坐标轴位置
        cb.Position([2, 4]) = [axPos(2), axPos(4)]; 
        
        % 复制坐标轴和Colorbar
        tmpFig = figure('Visible', 'off');
        new_handles = copyobj([ax, cb], tmpFig); % 向量形式复制关联对象
        new_ax = new_handles(1);
        new_cb = new_handles(2);
        
        % 在临时图窗中重调位置
        set(new_ax, 'Position', [0.1 0.1 0.7 0.8]);
        % 保持Colorbar与图像等高
        new_axPos = new_ax.Position;
        set(new_cb, 'Position', [new_axPos(1)+new_axPos(3)+0.03, new_axPos(2), 0.03, new_axPos(4)]); 

        switch k
            case 1
                if ~findFig(1)
                    saveas(tmpFig, fullfile(dir, figName{1})); 
                end
            case 2
                if ~findFig(2)
                    saveas(tmpFig, fullfile(dir, figName{2})); 
                end
            case 3
                saveas(tmpFig, fullfile(dir, figName{3})); 
        end
        close(tmpFig);
    end
    close(fig);
end