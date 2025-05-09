%% Plot a single result
clear;
K = 16;
r = 0.05;
method = "random_cmsbl";

file_name = sprintf("test/K=%d_r=%.2f_%s", K, r, method);
load(sprintf("test/%s.mat", file_name));

signal_origin = ln_to_db(phi * omega_real);
signal_recov = ln_to_db(phi * omega_est);

RSS_ideal = reshape(signal_origin, map.width, map.height);
RSS_real = reshape(signal_recov, map.width, map.height);

param = struct('mse', mse_db, 'sparsity', K, 'rate', r, 'method', method);

plot_recover_signal(map, RSS_ideal, RSS_real, param, file_name);

% Plot map information
function plot_recover_signal(map, RSS_ideal, RSS_real, param, file_name)
    % Params set
    Nx = map.width;
    Ny = map.height;

    % Plot part
    figure;
    subplot(1,3,1);
    imagesc(RSS_ideal);
    hold on; 
    axis([1 Nx 1 Ny]); 
    axis equal;
    axis tight;
    colorbar;
    title('Origin');
    xlabel(sprintf('K = %d, r = %.2f', param.sparsity, param.rate));
    ylabel('');
    ylabel(colorbar, 'dBm');
    hold off;
    
    subplot(1,3,2);
    imagesc(RSS_real);
    hold on; 
    axis([1 Nx 1 Ny]); 
    axis equal;
    axis tight;
    colorbar;
    title('Recover');
    xlabel(sprintf('mse = %.2f dB', param.mse));
    ylabel('');
    ylabel(colorbar, 'dBm');
    hold off;
    
    subplot(1,3,3);
    imagesc(RSS_ideal - RSS_real);
    hold on; 
    axis([1 Nx 1 Ny]); 
    axis equal;
    axis tight;
    colorbar;
    title('Diff');
    xlabel('');
    ylabel('');
    ylabel(colorbar, 'dBm');
    hold off;

    set(gcf, 'Position', [200, 100, 1200, 600]);
    set(gcf, 'Name', 'Origin vs Recover', 'NumberTitle', 'off');

    saveas(gcf,sprintf('%s.png', file_name));

end

% Convert linear to dB
function output = ln_to_db(input)
    output = 10 .* log10(input);
end