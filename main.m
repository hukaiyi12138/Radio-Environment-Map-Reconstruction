%% main
clear; close all; clc;

%% Set params
mapname         = "tinymap"    ; % {tinymap, largemap, hugemap}
K               = 4             ; % {4, 8, 12, 16}
sample_rate     = 0.15          ; % (sa) 0.01 ~ 0.15 / 0.05, 0.10, 0.15
select_rate     = 0.5           ; % (se) Select partial Rx
method_phi      = "idw"         ; % {idw, halrtc, kriging}
method_psi      = "random"      ; % {random, mmi}
method_recov    = "omp"         ; % {omp, sbl, csbl, msbl, cmsbl}
sigma2          = 0.05          ; % Noise power: σ^2

% Output directory
if ~exist("result", "dir")
    mkdir("result");
end

dir_out = sprintf("result/%s", mapname);
if ~exist(dir_out,"dir")
    mkdir(dir_out);
end

% Code start
fprintf('---- Code started ----\n');
method = sprintf('%s_%s_%s', method_phi, method_psi, method_recov); % entire method name
exp = sprintf('%s_K=%d_sa=%.2f_se=%.2f_si=%.2f_%s', mapname, K, sample_rate, select_rate, sigma2, method); % exper name
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
file_phi = sprintf('%s/%s_phi_%s_K=%d_se=%.2f.mat', dir_out, ...
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
mse = 10 * log10(norm(omega_real - omega_est) / norm(omega_real));
rmse = norm(rss_est_dbm(~isinf(rss_est_dbm)) - rss_full_dbm(~isinf(rss_full_dbm))) / sqrt(map.size);

fprintf('\n---- Evaluation: %s ----\n', exp);
fprintf('MSE = %.4f dB\n', mse);
fprintf('RMSE = %.4f dB\n', rmse);

%% Plot and save figure
plotRSSfig(map, omega_est, rss_full2D, rss_real2D, rss_est2D, ...
    sample_rate, mse, rmse, method, dir_out, exp);

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
function plotRSSfig(map, omega_est, rss_full2D, rss_real2D, rss_est2D, ...
    sample_rate, mse, rmse, method, dir_out, expName)
    % Find positions
    [buildX, buildY] = itc(map.build, map.Nx);
    
    % Set color parameters
    n = 512;
    h = linspace(0.7, 0, n)';
    s = linspace(0.8, 1, n)';
    v = linspace(0.9, 1, n)';
    cmap = hsv2rgb([h, s, v]);
    barMin = -75;
    barMax = 10;
    
    % Prepare figure and settings
    fig = figure('Visible', 'off', 'Position', [250, 100, 1200, 700]); % Suppress display
    titles = {'Ground Truth', 'RSS real', 'RSS est'};
    data = {rss_full2D, rss_real2D, rss_est2D};
    xlabels = {
        sprintf('Y axis\n\n%s', method), 
        sprintf('Y axis\n\nSelect Rate = %.0f%%\n\nSample Rate = %.2f', map.select_rate*100, sample_rate), 
        sprintf('Y axis\n\nMSE = %.4f dB\n\nRMSE = %.4f dB', mse, rmse)
    };
    
    % Generate and save each subplot separately
    for k = 1:3
        ax = subplot(1,3,k);
        imagesc(data{k});
        set(ax, 'YDir','normal');
        axis([1 map.Nx 1 map.Ny]);
        axis equal tight;
        hold(ax, 'on');
        colormap(ax, cmap);
        set(ax, 'CLim', [barMin barMax]);
        colorbar;
        ylabel(colorbar, 'dBm');
        title(titles{k});
        xlabel(xlabels{k}, 'Interpreter', 'none');
        ylabel('X axis');
        
        % Add markers
        dpi = get(0, 'ScreenPixelsPerInch');
        pos = ax.Position;
        cellW = pos(3) / map.Nx;
        cellH = pos(4) / map.Ny;
        markerPix = min(cellW, cellH);
        markerPt = markerPix * 72 / dpi;
        markerSize = markerPt^2 + 30;
        
        scatter(ax, buildY, buildX, markerSize, 's', ...
               'MarkerFaceColor','k', 'MarkerEdgeColor','none');
        
        % Save individual subplot
        tmpFig = figure('Visible', 'off');
        copyobj(ax, tmpFig);
        set(gca, 'Position', [0.1 0.1 0.8 0.8]); 
        saveas(tmpFig, fullfile(dir_out, sprintf('%s-%d.png', expName, k)));
        close(tmpFig);
    end
    
    close(fig); % Close main figure without saving
end