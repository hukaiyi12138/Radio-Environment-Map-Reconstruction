%% CSBL algorithm
% Cluster-based sbl algorithm
% y: Observation vector
% A: Sensing matrix
% t: Number of iterations
% sigma: Noise power (σ^2)
% mu: Recovered signal
function [mu, converge_point, mu_record] = csbl(y, A, t, sigma)
    % Initialize parameters
    [M, N] = size(A); % A: M x N
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    tol = 1e-3;

    % Calculate hyperparameters alpha and beta
    alpha = sigma^(-1) * ones(N, 1); % alpha: N x 1
    beta = sigma^(-1); % beta: scalar

    % Iteration process
    param1 = A' * A; % param1: N x N
    param2 = A' * y; % param2: N x 1

    Sigma = pinv(beta * param1 + diag(alpha));  % Sigma: N x N
    mu = beta * Sigma * param2; % N x 1

    mu_record = zeros(N, t);  % To store the value of mu in each iteration

    for iter = 1: t

        % Update alpha & beta
        alpha_new = (1 + 2 * a) ./ (mu.^2 + diag(Sigma) + 2 * b);
        beta_new = (M + 2 * c) / ((norm(y - A * mu))^2 + (beta^(-1)) * (N - alpha' * diag(Sigma)) + 2 * d);

        % Check for convergence
        converge_point = iter;
        if iter >= 2 && norm(mu - mu_record(:, iter - 1)) < tol
            converge_point = iter - 1;
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

    end
    
end
