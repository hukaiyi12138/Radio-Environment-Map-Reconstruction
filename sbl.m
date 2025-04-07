function [signal, alpha_cur, beta_cur] = sbl(y, A, iteration_time, sigma)
    % 初始化参数
    [M, N] = size(A); % A: M*N
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;

    % 计算超参数α和β
    alpha = sigma^(-1) * ones(N, 1); % alpha: N*1
    beta = sigma^(-1); % beta: 1*1

    % 变化曲线
    alpha_cur = 0;
    beta_cur = zeros(1, iteration_time);

    % 迭代过程
    param1 = A' * A; % param1: N*N
    param2 = A' * y; % param2: N*1
    Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N*N
    mu = beta * Sigma * param2; % N*1
    for iter = 1: iteration_time

        % 得到beta_cur
        beta_cur(iter) = beta;

        % 更新alpha & beta
        alpha_new = (1 + 2 * a) ./ (mu.^2 + diag(Sigma) + 2 * b);
        beta_new = (M + 2 * c)/(norm(y - A * mu)^2 + sigma * (1 - alpha' * diag(Sigma)) + 2 * d);
        alpha = alpha_new;
        beta = beta_new;

        % 计算协方差矩阵 Sigma 和均值 mu
        Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N*N
        mu = beta * Sigma * param2; % N*1

    end

    mu = max(mu, 0);

    % 恢复的信号
    signal = mu; % signal: N*1

end
