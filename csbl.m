%% CSBL algorithm
% Cluster-based sbl algorithm
% y: Observation vector
% A: Sensing matrix
% t: Number of iterations
% sigma: Noise power (σ^2)
% mu: Recovered signal
function [mu, iter] = csbl(y, A, t, sigma2)
    % Initialize parameters
    [M, N] = size(A); % A: M x N
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    tol_abs = 1e-4;
    res_old = -Inf;

    % Calculate hyperparameters alpha and beta
    alpha = sigma2^(-1) * ones(N, 1); % alpha: N x 1
    beta = sigma2^(-1); % beta: scalar

    % Iteration process
    param1 = A' * A; % param1: N x N
    param2 = A' * y; % param2: N x 1

    Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N x N
    mu = beta * Sigma * param2; % N x 1

    for iter = 1: t
        fprintf("SBL iter = %d\n", iter);

        % Update alpha & beta
        alpha_new = (1 + 2 * a) ./ (mu.^2 + diag(Sigma) + 2 * b);
        beta_new = (M + 2 * c) / ((norm(y - A * mu))^2 + (beta^(-1)) * (N - alpha' * diag(Sigma)) + 2 * d);

        % Check for convergence
        res_new = norm(y - A * mu)^2;
        if abs(res_new - res_old) <= tol_abs
            fprintf("SBL converge at %d-th iteration\n", iter);
            break;
        end

        alpha = alpha_new;
        beta = beta_new;

        % Adaptive dynamic threshold truncation (CSBL)
        threshold = mean(alpha.^(-1)) - std(alpha.^(-1));
        cond_alpha = alpha.^(-1) <= threshold;
        alpha(cond_alpha) = 0;  % Set small alpha values to zero
        alpha(alpha.^(-1) <= threshold) = 0;

        % Compute covariance matrix Sigma and mean mu
        Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N x N
        mu = beta * Sigma * param2; % N x 1

        % Update res_old
        res_old = res_new;
    end
end
