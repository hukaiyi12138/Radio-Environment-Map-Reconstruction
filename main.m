%% main
clear; close all; clc;

%% Set params
mapname         = "largemap"     ; % {tinymap, largemap, hugemap}
K               = 4             ; % {4, 8, 12, 16}；
sample_rate     = 0.15          ; % (sa) 0.01 ~ 0.15 / 0.05, 0.10, 0.15
select_rate     = 0.5           ; % (se) Select partial Rx
method_phi      = "idw"         ; % {idw, halrtc, kriging}
method_psi      = "random"      ; % {random, mmi}
method_recov    = "omp"         ; % {omp, sbl, csbl, msbl, cmsbl}
sigma           = 0.05          ; % Noise power: sigma = σ^2

% Output directory  
dir_out = sprintf("result/%s", mapname);
if ~exist(dir_out,"dir")
    mkdir(dir_out);
end

% Code start
fprintf('---- Code started ----\n');
method = sprintf('%s_%s_%s', method_phi, method_psi, method_recov); % entire method name
exp = sprintf('K=%d_sa=%.2f_se=%.2f_%s', K, sample_rate, select_rate, method); % exper name
fprintf('Map: %s\n', mapname);
fprintf('Method: %s\n', method);

% Load map file
file_map = sprintf('%s/%s.mat', dir_out, mapname);
if ~exist(file_map, "file")
    fprintf("Generating map: %s\n", mapname);
    map = generate_map2D(mapname);
    save(file_map, 'map');
else
    fprintf("Read map file: %s\n", file_map);
    load(file_map);
end

% Load map file with specific selected rate
file_map_se = sprintf('%s/%s_K=%d_se=%.2f.mat', dir_out, mapname, K, select_rate);
if ~exist(file_map_se, "file")
    fprintf("Selected rate = %.2f in map: %s\n",select_rate, mapname);

    % Select K Tx
    map.sparsity   = K;
    map.select_rate = select_rate;
    
    % Select partial Rx
    map.validPos = setdiff(1:map.size, map.build)'; % 除去建筑物的剩余位置
    map.validPosNum = numel(map.validPos);
    map.partPos = randsample(map.validPos, floor(select_rate * map.validPosNum)); % 对剩余位置进行随机采样
    map.interPos = setdiff(map.validPos, map.partPos); % 需要插值的位置
    
    % Construct ω_real
    map.selectedTxPos = [map.Tx(1:K).pos];
    map.omega_real = sparse(map.selectedTxPos, 1, [map.Tx(1:K).Pt_mW], map.size, 1);

    save(file_map_se, 'map');
else
    fprintf("Read map file: %s\n", file_map_se);
    load(file_map_se);
end

%% Generate sparse dictionary
file_phi = sprintf('%s/phi_%s_K=%d_se=%.2f.mat', dir_out, method_phi, K, select_rate);
if ~exist(file_phi, "file")
    fprintf("Generating φ by %s\n", method_phi);
    [phi, phi_rt, phi0] = generate_phi(method_phi, map);
    save(file_phi,'method_phi', 'phi', 'phi_rt', 'phi0');
else 
    fprintf("Read φ file: %s\n", file_phi);
    load(file_phi);
end

%% Generate measurement matrix
[psi] = generate_psi(method_psi, map, sample_rate, phi);

% Transmit process
omega_real = map.omega_real;
Phi = psi * phi; % Sensing matrix
y = Phi * omega_real; % Observation vector

%% Recover signal
[omega_est] = recover_signal(method_recov, y, Phi, sigma, K);

rss_full_vec = phi0 * map.omega_real;
rss_full_dbm = ln_to_db(rss_full_vec);
rss_full2D = reshape(rss_full_dbm, map.Nx, map.Ny);

rss_real_vec = phi_rt * map.omega_real;
rss_real_dbm = ln_to_db(rss_real_vec);
rss_real2D = reshape(rss_real_dbm, map.Nx, map.Ny);

rss_est_vec = phi * omega_est;
rss_est_dbm = ln_to_db(rss_est_vec);
rss_est2D = reshape(rss_est_dbm, map.Nx, map.Ny);

rss_full_dbm(rss_est_dbm == -Inf) = 0;
rss_real_dbm(rss_real_dbm == -Inf) = 0;
rss_est_dbm(rss_est_dbm == -Inf) = 0;

%% Evaluation
mse = 10 * log10(norm(omega_real - omega_est) / norm(omega_real));
rmse = norm(rss_est_dbm - rss_full_dbm) / sqrt(map.size);

fprintf('\n---- Evaluation: %s ----\n', exp);
fprintf('MSE = %.4f dB\n', mse);
fprintf('RMSE = %.4f dB\n', rmse);

%% Plot and save figure

% % plot phi_rt and phi
% figure("Visible","off");
% imagesc(phi_rt);
% axis([1 map.size 1 map.size]);
% axis equal;
% axis tight;
% colorbar;
% title(sprintf('φ^R^T - sparsity = %d', K));
% set(gcf, 'Position', [500, 200, 620, 500]);
% saveas(gcf, sprintf('%s/phiRT_K=%d.png', dir_out, K));
% 
% figure("Visible","off");
% imagesc(phi);
% axis([1 map.size 1 map.size]);
% axis equal;
% axis tight;
% colorbar;
% title(sprintf('%s', method_phi));
% set(gcf, 'Position', [500, 200, 620, 500]);
% saveas(gcf, sprintf('%s/phi_%s_K=%d.png', dir_out, method_phi, K));

% Tx and Rx information
% txX_real = [map.Tx.y];
% txY_real = [map.Tx.x];
[txY_real, txX_real] = itc(map.selectedTxPos, map.Nx);
[txY_est, txX_est] = itc(find(omega_est ~= 0), map.Nx);

% Contrast RSS figure
n = 512;
h = linspace(0.7, 0, n)';
s = linspace(0.8, 1, n)';
v = linspace(0.9, 1, n)';
cmap = hsv2rgb([h, s, v]);
barMin = -75;
barMax = 10;

figure;
subplot(1,3,1);
imagesc(rss_full2D);
set(gca,'YDir','normal');
axis([1 map.Nx 1 map.Ny]);
axis equal; axis tight; hold on;
% scatter(txX_real, txY_real, 30, 'white', 'filled');
colormap(gca, cmap);
set(gca, 'CLim', [barMin barMax]);
ylabel(colorbar, 'dBm');
title('Ground Truth');
xlabel(sprintf('Y axis\n\n%s', method), 'Interpreter','none');
ylabel('X axis');

subplot(1,3,2);
imagesc(rss_real2D);
set(gca,'YDir','normal');
axis([1 map.Nx 1 map.Ny]);
axis equal; axis tight; hold on;
% scatter(txX_real, txY_real, 30, 'white', 'filled');
colormap(gca, cmap);
set(gca, 'CLim', [barMin barMax]);
ylabel(colorbar, 'dBm');
title('RSS real');
xlabel(sprintf('Y axis\n\nSelect Rate = %.0f%%\n\nSample Rate = %.2f', map.select_rate * 100, sample_rate));
ylabel('X axis');

subplot(1,3,3);
imagesc(rss_est2D);
set(gca,'YDir','normal');
axis([1 map.Nx 1 map.Ny]);
axis equal; axis tight; hold on;
% scatter(txX_real, txY_real, 30, 'white', 'filled');
colormap(gca, cmap);
set(gca, 'CLim', [barMin barMax]);
ylabel(colorbar, 'dBm');
title('RSS est');
xlabel(sprintf('Y axis\n\nMSE = %.4f dB\n\nRMSE = %.4f dB', mse, rmse));
ylabel('X axis');

% subplot(1,4,4);
% imagesc(abs(rss_full2D - rss_est2D));
% set(gca,'YDir','normal');
% axis([1 map.Nx 1 map.Ny]);
% axis equal;
% axis tight;
% hold on;
% colorbar;
% ylabel(colorbar, 'dBm');
% title('Diff');
% xlabel('Y axis');
% ylabel('X axis');

set(gcf, 'Position', [250, 100, 1200, 700]);
drawnow;

% 标记建筑物
[buildX, buildY] = itc(map.build, map.Nx);      % 转成行列坐标
dpi = get(0,'ScreenPixelsPerInch');             % 屏幕 DPI，用于 px↔pt 转换
for k = 1:3
    ax = subplot(1,3,k);                        % 取出第 k 个子图
    ax.Units = 'pixels';
    pos = ax.Position;                          % [left bottom width height]（像素）
    
    % 每个格子的像素宽高
    cellW = pos(3) / map.Nx;
    cellH = pos(4) / map.Ny;
    markerPix = min(cellW, cellH);              % 选取较小的边长
    
    % 转成 pt，并计算面积（scatter 的 size 单位是 pt^2）
    markerPt = markerPix * 72 / dpi;
    markerSize = markerPt^2;
    
    hold(ax, 'on');
    scatter(ax, buildY, buildX, markerSize, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor','none');
    scatter(txX_real, txY_real, markerSize, 'white', 'filled');

end

saveas(gcf, sprintf('%s/rss_%s.png', dir_out, exp));

% save_rss_png_with_markers( ...
%     rss_full2D, cmap, barMin, barMax, ...
%     map.build, txX_real, txY_real, ...
%     map, fullfile(dir_out, sprintf('rss_%s_marked.png',exp)) );

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

% Draw RSS picture
function save_rss_png_with_markers( ...
        data,       ... % M×N 矩阵
        cmap,       ... % K×3 double [0,1]
        barMin, barMax, ...
        buildIdx,   ... % 建筑物在 data 展平后对应的线性索引向量
        txX, txY,   ... % 发射机在 data 矩阵中的行列坐标向量
        map,        ... % 包含 map.Nx（列数）信息，用于 itc 映射
        filename )  % 'xxx.png'
    
    %——1) 裁剪 & 归一化——
    data(data<barMin) = barMin;
    data(data>barMax) = barMax;
    normed = (data - barMin) / (barMax - barMin);
    
    %——2) 转成索引 & RGB8——
    K   = size(cmap,1);
    idx = uint16( floor(normed*(K-1)) ) + 1;
    RGB = ind2rgb(idx, cmap);           % M×N×3 double [0,1]
    RGB8 = uint8(RGB * 255);            % M×N×3 uint8
    
    %——3) 叠加“建筑物”黑色方块——
    [bX, bY] = itc(buildIdx, map.Nx);   % bX 行索引, bY 列索引
    ms = 30;  % 半边长像素，你可以根据需要调节
    for k = 1:numel(buildIdx)
        rows = max(1,bX(k)-ms) : min(size(RGB8,1), bX(k)+ms);
        cols = max(1,bY(k)-ms) : min(size(RGB8,2), bY(k)+ms);
        RGB8(rows,cols,:) = 0;  % 全黑
    end
    
    %——4) 叠加“发射机”白色点——
    for k = 1:numel(txX)
        x = round(txX(k)); 
        y = round(txY(k));
        if x>=1 && x<=size(RGB8,1) && y>=1 && y<=size(RGB8,2)
            RGB8(x,y,:) = 255;    % 单像素白点
        end
    end
    
    %——5) 写出带标注的彩色 PNG——
    imwrite(RGB8, filename);
end

