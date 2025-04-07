%% Plot basic information
function plot_basic(map_index)
    % Load data
    fileName = sprintf('output/map_%d/map%d.mat', map_index, map_index);
    load(fileName);

    %%% Figure 1: RSS_entire 和 RSS_ideal 对比图 %%%
    figure('Visible', 'off');
    subplot(1,2,1);  % 第一张图放在窗口的左侧
    imagesc(RSS_entire);
    hold on; 
    set(gcf, 'Name', 'Figure 1: RSS Comparison', 'NumberTitle', 'off');
    axis([1 Nxx 1 Nyy]); 
    axis equal;  % 设置为正方形像素
    axis tight;
    scatter(Tx.X, Tx.Y, 20, 'r', 'filled');  % 绘制发射机位置
    rectangle('Position', [Nx1, Ny1, Nx2-Nx1, Ny2-Ny1], 'EdgeColor', 'w', 'LineWidth', 2); % 框出ROI
    colorbar; % 添加颜色条以显示数值范围
    title('RSS entire');
    xlabel('X axis');
    ylabel('Y axis');
    ylabel(colorbar, 'dBm');
    hold off;
    
    subplot(1,2,2);  % 第二张图放在窗口的右侧
    imagesc(RSS_ideal);
    hold on; 
    axis([1 Nx 1 Ny]); 
    axis equal;  % 设置为正方形像素
    axis tight;
    scatter(Tx_roi.X, Tx_roi.Y, 20, 'r', 'filled');  % 绘制发射机位置
    colorbar; % 添加颜色条以显示数值范围
    title('RSS ideal');
    xlabel('X axis');
    ylabel('Y axis');
    ylabel(colorbar, 'dBm');
    hold off;
    
    set(gcf, 'Position', [100, 100, 1200, 600]);  % 设置窗口大小和位置
    saveas(gcf, sprintf('output/map_%d/map_%d-RSS_entire.png', map_index ,map_index));
    
    %%% Figure 2: phi_rt 和 phi 对比图 %%%
    figure('Visible', 'off');
    subplot(1,2,1);  % 第一张图放在窗口的左侧
    imagesc(phi_rt);
    set(gcf, 'Name', 'Figure 2: Sparse Dictionary Design', 'NumberTitle', 'off');
    axis([1 N 1 N]);
    axis equal;  % 设置为正方形像素
    axis tight;
    colorbar;  % 添加颜色条以显示数值范围
    title('phi RT');
    xlabel('X axis');
    ylabel('Y axis');
    
    subplot(1,2,2);  % 第二张图放在窗口的右侧
    imagesc(phi);
    axis([1 N 1 N]);
    axis equal;  % 设置为正方形像素
    axis tight;
    colorbar;  % 添加颜色条以显示数值范围
    title('phi');
    xlabel('X axis');
    ylabel('Y axis');
    
    set(gcf, 'Position', [100, 100, 1200, 600]);  % 设置窗口大小和位置
    saveas(gcf, sprintf('output/map_%d/map_%d-Sparse_Dictionary.png', map_index, map_index));

end