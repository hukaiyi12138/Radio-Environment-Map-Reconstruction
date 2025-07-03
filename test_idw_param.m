%% Test: IDW evaluation according to different parameter p
clear; close all; clc;

%% Set params
mapname         = "hugemap"     ; % {tinymap, largemap, hugemap}
K               = 4             ; % {4, 8, 12, 16}；
select_rate     = 0.5           ; % (se) Select partial Rx
method_phi      = "idw"         ; % {idw, halrtc, kriging}

% Output directory  
exp = sprintf('%s_%s_K=%d_se=%.2f', method_phi, mapname, K, select_rate); % exper name
dir_in = sprintf("result/%s", mapname);
if ~exist(dir_in,"dir")
    error("directory %s does not exist.\n", dir_in);
end
dir_out = sprintf("test/%s", exp);
if ~exist(dir_out,"dir")
    mkdir(dir_out);
end

% Load map file with specific selected rate
file_map_se = sprintf('%s/%s_K=%d_se=%.2f.mat', dir_in, mapname, K, select_rate);
if ~exist(file_map_se, "file")
    error("Map info does not exist.\n");
else
    fprintf("Read map file: %s\n", file_map_se);
    load(file_map_se);
end

%% Generate sparse dictionary
diary 'hugemap.txt';

% 设置一系列 p 值
p_list = 5:25;
rmse_list = zeros(size(p_list));

% 对每个 p 计算 rmse
for i = 1:numel(p_list)
    p = p_list(i);

    % 1) 重新生成稀疏字典 φ
    [phi_rt, phi0] = generate_phirt(map);
    phi = phi_idw(phi_rt, map, p);

    % 2) 计算全量和估计的 RSS（dBm）
    rss_full_vec = phi0 * map.omega_real;
    rss_full_dbm = ln_to_db(rss_full_vec);
    rss_est_vec  = phi  * map.omega_real;
    rss_est_dbm  = ln_to_db(rss_est_vec);

    % 3) 重塑为 2D 矩阵
    rss_full2D = reshape(rss_full_dbm, map.Nx, map.Ny);
    rss_est2D  = reshape(rss_est_dbm,  map.Nx, map.Ny);
    rss_full_dbm2 = rss_full_dbm;
    rss_est_dbm2 = rss_est_dbm;

    % 计算 rmse
    rss_full_dbm2(rss_est_dbm == -Inf) = 0;
    rss_est_dbm2(rss_est_dbm == -Inf) = 0;
    rmse = norm(rss_est_dbm2 - rss_full_dbm2) / sqrt(map.size);
    rmse_list(i) = rmse;
    fprintf('p = %d, RMSE = %.4f\n', p, rmse_list(i));

    % 4) 绘图
    figure('Visible','off');
    n = 512;
    h = linspace(0.7,0,n)';
    s = linspace(0.8,1,n)';
    v = linspace(0.9,1,n)';
    cmap = hsv2rgb([h s v]);
    barMin = -75; barMax = 10;

    subplot(1,2,1);
    imagesc(rss_full2D);
    set(gca,'YDir','normal');
    axis equal tight; hold on;
    colormap(gca, cmap);
    c = colorbar; ylabel(c,'dBm');
    clim([barMin barMax]);
    title('Ground Truth');
    xlabel(sprintf('Y axis\n\np = %d', p));
    ylabel('X axis');

    subplot(1,2,2);
    imagesc(rss_est2D);
    set(gca,'YDir','normal');
    axis equal tight; hold on;
    colormap(gca, cmap);
    c = colorbar; ylabel(c,'dBm');
    clim([barMin barMax]);
    title('Interpolated');
    xlabel(sprintf('Y axis\n\nrmse = %.4f dB', rmse));
    ylabel('X axis');

    set(gcf, 'Position', [250, 100, 900, 700]);

    % 5) 保存
    fname = sprintf('idw_%s_p=%02d_K=%d_se=%.2f.png', mapname, p, K, select_rate);
    saveas(gcf, fullfile(dir_out, fname));
    close(gcf);

    fprintf('已保存：%s\n', fullfile(dir_out, fname));

    % 找出最优 p 值
    if i == 1
        best.p = p;
        best.rmse = rmse;
        best.phi0 = phi0;
        best.phi_rt = phi_rt;
        best.phi = phi;
    end
    if rmse < best.rmse
        best.p = p;
        best.rmse = rmse;
        best.phi0 = phi0;
        best.phi_rt = phi_rt;
        best.phi = phi;
    end
end

%% Plot RMSE curve
figure;
plot(p_list, rmse_list, '-o', 'LineWidth', 1.5);
hold on;
grid on;

% 找到最小 RMSE 及对应 p
[min_rmse, idx_min] = min(rmse_list);
p_opt = p_list(idx_min);
fprintf('最小 RMSE = %.4f， 对应 p = %d\n', min_rmse, p_opt);

% 在图中标出最小点
plot(p_opt, min_rmse, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% 添加文字注释
text(p_opt, min_rmse, sprintf('  p=%d\n  RMSE=%.2f', p_opt, min_rmse), ...
     'VerticalAlignment', 'bottom', 'Color', 'red', 'FontWeight', 'bold');

xlabel('Param p');
ylabel('RMSE (dB)');
title(sprintf('IDW 方法 RMSE 随 p 变化 (%s K=%d se=%.2f)', mapname, K, select_rate));

hold off;

% 保存图像
saveas(gcf, fullfile(dir_out, sprintf('idw_%s_K=%d_se=%.2f.png', mapname, K, select_rate)));
fprintf('已保存带标注的 RMSE 曲线到 %s\n', fullfile(dir_out, sprintf('%s.png', mapname)));

% 保存最优 p 值时的字典
phi0 = best.phi0;
phi_rt = best.phi_rt;
phi = best.phi;
file_phi = sprintf('%s/phi_%s_K=%d_se=%.2f.mat', dir_in, method_phi, K, select_rate);
save(file_phi,'method_phi', 'phi', 'phi_rt', 'phi0');
fprintf('--> Best p = %d <--\n', p);

diary off;

%% Core function
function [phi_rt, phi0] = generate_phirt(map)
    N    = map.size;
    
    phi0 = sparse(N, N);
    phi0(:, map.selectedTxPos) = map.phi(:, map.selectedTxPos);

    phi_rt = phi0;
    phi_rt(map.interPos, :) = 0;
end

%% Additional function
% Convert linear to dB
function db = ln_to_db(ln)
    db = 10.*log10(ln);
end
