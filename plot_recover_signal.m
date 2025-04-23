%% Plot recovered signal
function plot_recover_signal(map, )

    figure;
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
    xlabel(sprintf('MSE = %.4f dB', mse_random_omp));
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
    xlabel(sprintf('MSE = %.4f dB', mse_random_sbl));
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
    xlabel(sprintf('MSE = %.4f dB', mse_mmi_sbl));
    ylabel(''); 
    ylabel(colorbar, 'dBm');
    hold off;
    
    set(gcf, 'Position', [350, 150, 900, 600]);  % 设置窗口大小和位置
    set(gcf, 'Name', sprintf('Recovered Map - Sampling Rate=%.2f', rate), 'NumberTitle', 'off');

end
