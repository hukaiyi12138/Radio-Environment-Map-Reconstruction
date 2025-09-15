%% main
clear; close all; clc;

%% Set params
mapname         = "tinymap"    ; % {tinymap, largemap, hugemap, bupt}
K               = 4             ; % {4, 8, 12, 16}
sample_rate     = 0.15          ; % (sa) 0.01 ~ 0.15 / 0.05, 0.10, 0.15
select_rate     = 0.5           ; % (se) Select partial Rx
sigma2          = 0.05          ; % Noise power: σ^2
method_phi      = "idw"         ; % {idw, halrtc, kriging}
method_psi      = "random"      ; % {random, mmi}
method_recov    = "omp"         ; % {omp, sbl, csbl, msbl, cmsbl}
figPlot         = true          ; % plot figures or not

% Output directory
if ~exist("result", "dir")
    mkdir("result");
end

dir = sprintf("result/%s", mapname); % map path
if ~exist(dir,"dir")
    mkdir(dir);
end

dir_fig = fullfile(dir, "figures"); % figures path
if ~exist(dir_fig, "dir")
    mkdir(dir_fig); 
end

% Code start
fprintf('---- Code started ----\n');
method = sprintf('%s_%s_%s', method_phi, method_psi, method_recov); % entire method name
exp = sprintf('%s_%s_K=%d_sa=%.2f_se=%.2f_si=%.2f', mapname, method, K, sample_rate, select_rate, sigma2);
fprintf('Map: %s\n', mapname);
fprintf('Method: %s\n', method);

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
file_psi = sprintf('%s/%s_psi_%s_%s_K=%d_sa=%.2f_se=%.2f.mat', dir, ...
        mapname, method_psi, method_phi, K, sample_rate, select_rate);
if ~exist(file_psi, "file") 
    fprintf("Generating φ by %s\n", method_psi);
    [psi] = generate_psi(method_psi, map, sample_rate, phi);
    if method_psi == "mmi"
        save(file_psi,'method_psi', 'psi');
    end
else 
    fprintf("Read ψ file: %s\n", file_psi);
    load(file_psi);
end

% Transmit process
omega_real = map.omega_real;
Phi = psi * phi; % Sensing matrix
noise = sqrt(sigma2) * randn(size(Phi, 1), 1); % Gaussian noise
y = Phi * omega_real + noise; % Observation vector

%% Recover signal
result_name = sprintf('%s/%s.mat', dir, exp);
if ~exist(result_name, 'file')
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
else
    fprintf("Read exp file: %s\n", result_name);
    load(result_name);
end

%% Evaluation
if ~exist(result_name, 'file')
    valid.mask = isfinite(rss_est_dbm) & isfinite(rss_full_dbm); % ignore Inf value
    valid.num = numel(valid.mask);
    
    % MSE - sparse signal recovery
    mse = 10 * log10(norm(omega_real - omega_est) / norm(omega_real));
    % RMSE - RSS recovery
    rmse = norm(rss_est_dbm(valid.mask) - rss_full_dbm(valid.mask)) / sqrt(valid.num);
    % MAE - dictionary
    mae = norm(rss_est_dbm(valid.mask) - rss_full_dbm(valid.mask), 1) / valid.num;

    clear("valid");
end

fprintf('\n---- Evaluation: %s ----\n', exp);
fprintf('MSE = %.4f dB\n', mse);
fprintf('RMSE = %.4f dB\n', rmse);
fprintf('MAE = %.4f dB\n', mae);

%% Plot and save figure
% Plot figures
figName = {
    sprintf('%s-Ground_Truth.png', map.name), ...
    sprintf('%s_K=%d_se=%.2f.png', map.name, K, select_rate), ...
    sprintf('%s.png', exp)
};

% draw figs
if figPlot
    plotRSSfig(map, method, rss_full2D, rss_real2D, rss_est2D, figName, dir_fig);
end

% Save result
if ~exist(result_name, 'file')
    save(result_name);

    % data repository
    params = struct( ...
        'map',              mapname, ...
        'exper',            string(exp), ...
        'method',           string(method), ...
        'method_phi',       method_phi, ...
        'method_psi',       method_psi, ...
        'method_recovery',  method_recov, ...
        'sparsity',         K, ...
        'sample_rate',      sample_rate, ...
        'select_rate',      select_rate, ...
        'noise_power',      sigma2 ...
    );

    metrics = struct('MSE', mse, 'RMSE', rmse, 'MAE', mae);

    logfile = fullfile(dir, 'ExperLog.csv');
    log_experiment(logfile, params, metrics);
end

fprintf('\n---- Code finished ----\n');

%% Additional function
% Convert dB to linear
% function ln = db_to_ln(db)
%     ln = 10.^(db./10);
% end

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
function plotRSSfig(map, method, rss_full2D, rss_real2D, rss_est2D, figName, dir)
    [buildX, buildY] = itc(map.build, map.Nx);

    % 颜色和范围
    n = 512;
    cmap = hsv2rgb([linspace(0.7,0,n)', linspace(0.8,1,n)', linspace(0.9,1,n)']);
    barMin = -75; barMax = 10;

    % 跳过：三张都存在就不画
    paths = cellfun(@(f) fullfile(dir, f), figName, 'UniformOutput', false);
    if all(cellfun(@isfile, paths)), return; end

    fig = figure('Visible','off','Position', [250,100,2000,1000],'Color','w');
    titles = {'Ground Truth', 'Sparse Sampling', sprintf('%s', method)};
    data   = {rss_full2D, rss_real2D, rss_est2D};

    for k = 1:3
        D = data{k};
        D(isinf(D) & D < 0) = NaN;

        ax = subplot(1,3,k);
        set(ax,'Color','w');
        imagesc(ax, D, 'AlphaData', ~isnan(D));
        set(ax,'YDir','normal'); axis([1 map.Nx 1 map.Ny]); axis equal tight;
        colormap(ax, cmap); set(ax,'CLim',[barMin barMax]);

        cb = colorbar(ax,'Location','eastoutside'); ylabel(cb,'dBm');
        title(titles{k}, 'Interpreter','none');
        xlabel('Y axis'); ylabel('X axis');

        % 快速建筑覆盖
        overlay_buildings(ax, map, buildX, buildY);

        % 导出该子图
        out_path = fullfile(dir, figName{k});
        need_save = (k < 3 && ~isfile(out_path)) || k == 3;
        if need_save
            exportgraphics(ax, out_path, 'Resolution', 180);
        end
    end
    close(fig);
end

% Plot building positions
function overlay_buildings(ax, map, buildX, buildY)
    B = false(map.Nx, map.Ny);
    B(sub2ind([map.Nx, map.Ny], buildX, buildY)) = true;
    C = zeros(map.Nx, map.Ny, 3);  % 黑色
    hold(ax,'on');
    hi = image(ax, C);
    set(hi, 'AlphaData', B);
    uistack(hi,'top');
end

% Save exper data
function log_experiment(logfile, params, metrics)
    % % 拼一行
    % row = struct();
    % row.timestamp = string(datetime('now','TimeZone','local','Format','yyyy-MM-dd HH:mm:ss'));
    % 参数
    pf = fieldnames(params);
    for i = 1:numel(pf); row.(pf{i}) = params.(pf{i}); end
    % 指标
    mf = fieldnames(metrics);
    for i = 1:numel(mf); row.(mf{i}) = metrics.(mf{i}); end

    Trow = struct2table(row); % 标准化为 table

    if isfile(logfile)
        Texist = readtable(logfile);

        % 类型/列名对齐（以新行 Trow 为模板）
        Texist = harmonize_types_and_vars(Texist, Trow);
        Trow   = harmonize_types_and_vars(Trow, Texist); % 双向一次，确保顺序一致

        % 以 exp 为主键 upsert（没有 exp 列则纯追加）
        if ismember('exp', Texist.Properties.VariableNames)
            exp_exist = string(Texist.exp);
            exp_new   = string(Trow.exp);
            idx = strcmp(exp_exist, exp_new);
        else
            idx = false(height(Texist),1);
        end

        if any(idx)
            Texist(idx, :) = Trow(:, Texist.Properties.VariableNames);
        else
            Texist = [Texist; Trow(:, Texist.Properties.VariableNames)];
        end

        writetable(Texist, logfile);
    else
        % 首次写入
        writetable(Trow, logfile);
    end
end

function T = harmonize_types_and_vars(T, Tref)
    % 确保 T 至少拥有 Tref 的所有列，且类型与 Tref 对齐（double 或 string）
    refVars = Tref.Properties.VariableNames;

    % 先补列
    for k = 1:numel(refVars)
        vn = refVars{k};
        if ~ismember(vn, T.Properties.VariableNames)
            % 依据 Tref 的类型选择默认值
            cls = class(Tref.(vn));
            if ismember(cls, {'double','single'})
                T.(vn) = nan(height(T),1);
            elseif ismember(cls, {'string','char','cell'})
                T.(vn) = strings(height(T),1);
            else
                % 其他类型一律当字符串存
                T.(vn) = strings(height(T),1);
            end
        end
    end

    % 再做类型对齐
    for k = 1:numel(refVars)
        vn  = refVars{k};
        cls_ref = class(Tref.(vn));
        cls_cur = class(T.(vn));
        if ~strcmp(cls_ref, cls_cur)
            try
                if ismember(cls_ref, {'double','single'})
                    % 目标是数值
                    if iscell(T.(vn)); T.(vn) = string(T.(vn)); end
                    if isstring(T.(vn)) || ischar(T.(vn))
                        T.(vn) = double(str2double(string(T.(vn))));
                    elseif islogical(T.(vn))
                        T.(vn) = double(T.(vn));
                    else
                        T.(vn) = double(T.(vn));
                    end
                else
                    T.(vn) = string(T.(vn));
                end
            catch
                % 若转换失败，回退为字符串
                T.(vn) = string(T.(vn));
            end
        end
    end

    % 变量顺序按参考表
    T = T(:, refVars);
end
