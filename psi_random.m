function psi = psi_random(M, N, sample_rate)
    % M, N: psi 矩阵的尺寸
    % 生成一个全零矩阵
    psi = zeros(M, N);
    
    % 确保 target_ones 不大于矩阵的元素总数
    samples = round(sample_rate * N);
    if samples > M * N
        error('target_ones cannot be greater than the total number of elements in the matrix.');
    end
    
    % 对每一行进行处理
    for i = 1:M
        ones_positions = randperm(N, samples);  % 随机选择位置
        psi(i, ones_positions) = 1;  % 将对应位置设置为 1
    end
end
