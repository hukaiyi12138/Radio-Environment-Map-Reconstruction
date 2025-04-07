%% Plot MSE vs. Sampling_rate
function plot_mse(map)
    load(sprintf('output/map_%d/result%d.mat', map, map));

    % Figure 4: 绘制每种算法的 MSE 与 sample_rate 的关系图
    mse_random_omp = arrayfun(@(x) x.mse, result_random_omp);  % Extract mse for Random OMP
    mse_random_sbl = arrayfun(@(x) x.mse, result_random_sbl);  % Extract mse for Random SBL
    mse_mmi_sbl = arrayfun(@(x) x.mse, result_mmi_sbl);  % Extract mse for Random SBL
    
    % Plot the MSE vs. Sampling rate
    figure;
    hold on;
    plot(sample_rate_values, mse_random_omp, '-o', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', 'Random OMP');
    plot(sample_rate_values, mse_random_sbl, '-s', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Random SBL');
    plot(sample_rate_values, mse_mmi_sbl, '-d', 'MarkerFaceColor', 'm', 'MarkerSize', 4, 'DisplayName', 'MMI SBL');
    hold off;
    
    xlabel('Sampling rate');
    ylabel('Sparse signal recovery error');
    title('The MSE of sparse signal recovery performance comparisons');
    legend('show');
    grid on;
    xlim([0.01, 0.15]); % 设置x轴范围从0.01到0.15
    
    set(gcf, 'Position', [350, 150, 900, 600]);  % 设置窗口大小和位置
    saveas(gcf, sprintf('output/map_%d/map_%d-MSE_vs_SamplingRate.png', map, map));

end