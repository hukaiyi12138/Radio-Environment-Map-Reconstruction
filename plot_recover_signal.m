%% Plot recovered signal
function plot_recover_signal(map, targaet_index)
    % Load data
    load(sprintf('output/map_%d/map%d.mat', map, map), 'phi', 'Nx', 'Ny', 'RSS_ideal');
    load(sprintf('output/map_%d/result%d.mat', map, map));
    rate = sample_rate_values(targaet_index);
    
    % Get mse values from the results
    mse_random_omp = result_random_omp(targaet_index).mse;  % MSE for Random-OMP
    mse_random_sbl = result_random_sbl(targaet_index).mse;  % MSE for Random-SBL
    mse_mmi_sbl = result_mmi_sbl(targaet_index).mse;  % MSE for MMI-SBL
    
    % Calculate RSS_real
    RSS_real1 = trans_rss(result_random_omp(targaet_index), phi, Nx, Ny);
    RSS_real2 = trans_rss(result_random_sbl(targaet_index), phi, Nx, Ny);
    RSS_real3 = trans_rss(result_mmi_sbl(targaet_index), phi, Nx, Ny);
    
    %%% Figure 3: RSS_ideal 和 RSS_real 对比图 %%%
    figure('Visible', 'off');
    tiledlayout(2,2);  % 使用 tiledlayout 替代 subplot
    
    % First plot: RSS_ideal
    nexttile;
    imagesc(RSS_ideal);
    hold on;
    axis([1 Nx 1 Ny]);
    axis equal;
    axis tight;
    colorbar;
    title('Original');
    xlabel('');
    ylabel('');
    ylabel(colorbar, 'dBm');
    hold off;
    
    % Second plot: RSS_real1 (Random-OMP)
    nexttile;
    imagesc(RSS_real1);
    hold on;
    axis([1 Nx 1 Ny]);
    axis equal;
    axis tight;
    colorbar;
    title('Random-OMP');
    xlabel(sprintf('MSE = %.2f dB', mse_random_omp));
    ylabel(''); 
    ylabel(colorbar, 'dBm');
    hold off;
    
    % Third plot: RSS_real2 (Random-SBL)
    nexttile;
    imagesc(RSS_real2);
    hold on;
    axis([1 Nx 1 Ny]);
    axis equal;
    axis tight;
    colorbar;
    title('Random-SBL');
    xlabel(sprintf('MSE = %e dB', mse_random_sbl));
    ylabel(''); 
    ylabel(colorbar, 'dBm');
    hold off;
    
    % Fourth plot: RSS_real3 (MMI-SBL)
    nexttile;
    imagesc(RSS_real3);
    hold on;
    axis([1 Nx 1 Ny]);
    axis equal;
    axis tight;
    colorbar;
    title('MMI-SBL');
    xlabel(sprintf('MSE = %e dB', mse_mmi_sbl));
    ylabel(''); 
    ylabel(colorbar, 'dBm');
    hold off;
    
    set(gcf, 'Position', [350, 150, 900, 600]);  % 设置窗口大小和位置
    set(gcf, 'Name', sprintf('Recovered Map - Sampling Rate=%.2f', rate), 'NumberTitle', 'off');
%     saveas(gcf, sprintf('output/map_%d/Recovered_Map-rate=%.2f.png', map, rate));

end

function RSS_real = trans_rss(result, phi, Nx, Ny)
    omega_est = result.omega_est;
    x_recov = phi * omega_est;
    x_recov = 10 * log10(x_recov);  % 转换为 dBm
    RSS_real = reshape(x_recov, Nx, Ny);
end
