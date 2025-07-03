%% SBL algorithm
% y: Observation vector
% A: Sensing matrix
% t: Number of iterations
% sigma: Noise power (σ^2)
% mu: Recovered signal
function [mu, iter] = sbl(y, A, t, sigma)
    % Initialize parameters
    [M, N] = size(A); % A: M x N
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    tol_abs = 1e-6;
    tol_rel = 1e-4; 

    % Calculate hyperparameters alpha and beta
    alpha = sigma^(-1) * ones(N, 1); % alpha: N x 1
    beta = sigma^(-1); % beta: scalar

    % Iteration process
    param1 = A' * A; % param1: N x N
    param2 = A' * y; % param2: N x 1

    Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N x N
    mu = beta * Sigma * param2; % N x 1

    res_prev = norm(y)^2;
    alpha_thre = 1e-4;
    for iter = 1: t

        % Update alpha & beta
        alpha_new = (1 + 2 * a) ./ (mu.^2 + diag(Sigma) + 2 * b);
        beta_new = (M + 2 * c) / ((norm(y - A * mu))^2 + (beta^(-1)) * (N - alpha' * diag(Sigma)) + 2 * d);

        % Check for convergence
        res = norm(y - A*mu)^2;
        if res < tol_abs || abs(res - res_prev)/res_prev < tol_rel
            break;
        end
        res_prev = res;

        alpha = alpha_new;
        beta = beta_new;

        % alpha threshold
        cond_alpha = alpha.^(-1) <= alpha_thre;
        alpha(cond_alpha) = 0;  % Set small alpha values to zero

        % Compute covariance matrix Sigma and mean mu
        Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N x N
        mu = beta * Sigma * param2; % N x 1
%         mu(abs(mu) < tol_abs) = 0;

    end
end
