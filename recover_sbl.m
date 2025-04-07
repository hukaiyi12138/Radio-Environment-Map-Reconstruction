%% Recover Signal by SBL
function [omega_est, mse, beta_cur] = recover_sbl(phi, psi, omega_real, noise, iteration_time, sigma)
    % 计算感知矩阵 A 和观测向量 y
    A = psi * phi;  % 感知矩阵 A
    y = A * omega_real + noise;  % 观测向量 y
    
    % 使用 SBL 算法恢复信号
    [omega_est, ~, beta_cur] = sbl(y, A, iteration_time, sigma);  % 调用 SBL 函数并获得恢复结果

    % MMD: 基于最大最小距离的聚类算法
    
    
    % 计算恢复信号的 MSE: 均方误差 (dB)
    mse = norm(omega_real - omega_est) / norm(omega_real);
    mse = 10 * log10(mse);  % 转换为dB
    
end
