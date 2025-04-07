%% Recover Signal by OMP
function [omega_est, mse, rn_norm, converge_point] = recover_omp(phi, psi, omega_real, noise, iteration_time)
    % 计算感知矩阵 A 和观测向量 y
    A = psi * phi;  % 感知矩阵 A
    y = A * omega_real + noise;  % 观测向量 y
    
    % 使用 OMP 算法恢复信号
    [omega_est, rn_norm, converge_point] = omp(y, A, iteration_time);

    % 计算 MSE: 均方误差 (dB)
    mse = norm(omega_real - omega_est) / norm(omega_real);
    mse = 10 * log10(mse);
    
end

