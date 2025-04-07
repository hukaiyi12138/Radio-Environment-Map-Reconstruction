function [psi, psi_opt] = generate_psi(sample_rate_values, M, N, phi)
    % 初始化 psi 为 cell 数组
    psi = cell(1, length(sample_rate_values));  % 预分配，存储不同采样率下的结果
    psi_opt = cell(1, length(sample_rate_values));  % 预分配用于存储 MMI 结果

    for i = 1:length(sample_rate_values)
        % 使用 cell 数组时要使用 `{}` 来访问
        psi{i} = psi_random(M, N, sample_rate_values(i));  
        fprintf('Finish psi_random: %d/%d\n', i, length(sample_rate_values));
        psi_opt{i} = psi_mmi(M, N, sample_rate_values(i), phi);  % 计算 MMI 采样矩阵
        fprintf('Finish psi_mmi: %d/%d\n', i, length(sample_rate_values));
    end
end
