function [mu, iter] = sbl(y, A, t, sigma2, algo)
% SBL  Sparse Bayesian Learning (SBL) / Cluster-based SBL (CSBL)
%
%   [mu, iter] = SBL(y, A, t, sigma2, algo)
%
%   Implements the Sparse Bayesian Learning algorithm and optionally the
%   cluster-based variant (CSBL). Uses a Sigma operator based on the
%   Woodbury identity to avoid explicitly storing the large covariance
%   matrix.
%
%   INPUT:
%     y      - (M x 1) observation vector
%     A      - (M x N) sensing matrix
%     t      - maximum number of iterations
%     sigma2 - noise variance (σ^2)
%     algo   - 'true' or true: use CSBL
%              otherwise: use standard SBL
%
%   OUTPUT:
%     mu     - (N x 1, sparse) recovered signal
%     iter   - number of iterations performed
%
%   NOTES:
%     - The covariance matrix Sigma = (β A'A + diag(α))^(-1) is not
%       explicitly formed (size N x N). Instead, a Sigma operator is
%       provided with apply/diag/cols interfaces.
%     - Hyperparameters α and β are updated using standard SBL/CSBL rules.
%
%   See also: chol, pcg

    if nargin < 5, algo = false; end
    if ischar(algo) || isstring(algo), algo = strcmpi(string(algo), 'true'); end

    [M, N] = size(A);

    % ---------------------------------------------------------------------
    % Hyperprior parameters
    % ---------------------------------------------------------------------
    a = 1e-6; b = 1e-6; c = 1e-6; d = 1e-6;
    tol_abs = 1e-4;          % absolute tolerance for convergence
    res_old = -Inf;
    alpha_thre = 1e-2;       % pruning threshold for SBL
    eps_alpha  = 1e-12;      % lower bound for alpha (to avoid zeros)

    % ---------------------------------------------------------------------
    % Initialization
    % ---------------------------------------------------------------------
    alpha = (1/sigma2) * ones(N,1);   % (N x 1) ARD precision
    beta  = 1/sigma2;                 % noise precision

    ATy = A.' * y;

    % Construct Sigma operator and initial estimates
    Sigma = makeSigmaOp(A, max(alpha, eps_alpha), beta);
    mu    = Sigma.apply(beta * ATy);      % mu = beta * Sigma * ATy
    dSig  = Sigma.diag();                 % diag(Sigma)

    % ---------------------------------------------------------------------
    % Iterative updates
    % ---------------------------------------------------------------------
    for iter = 1:t
        fprintf("SBL iter = %d\n", iter);

        % --- Hyperparameter updates ---
        alpha_new = (1 + 2*a) ./ (mu.^2 + dSig + 2*b);

        resid2    = norm(y - A*mu)^2;
        gamma_eff = sum(1 - alpha_new .* dSig);
        beta_new  = (M + 2*c) / ( resid2 + (1/beta)*gamma_eff + 2*d );

        % --- Convergence check ---
        res_new = resid2;
        if abs(res_new - res_old) <= tol_abs
            fprintf("SBL converged at iteration %d\n", iter);
            break;
        end
        res_old = res_new;

        % --- Sparsification step ---
        if algo
            % CSBL: adaptive threshold (mean - std of 1/alpha)
            inva = 1 ./ alpha_new;
            threshold = mean(inva) - std(inva);
            cond = inva <= threshold;
        else
            % Standard SBL: fixed threshold
            cond = (1 ./ alpha_new) <= alpha_thre;
        end
        alpha_new(cond) = max(alpha_new(cond), eps_alpha);

        % --- Apply updates ---
        alpha = alpha_new;
        beta  = beta_new;

        % --- Reconstruct Sigma and update estimates ---
        Sigma = makeSigmaOp(A, max(alpha, eps_alpha), beta);
        mu    = Sigma.apply(beta * ATy);
        dSig  = Sigma.diag();
    end

    mu = sparse(mu);
end

% ======================================================================
%                   Sigma Operator (Woodbury identity)
% ======================================================================
function Sigma = makeSigmaOp(A, alpha, beta)
% makeSigmaOp Construct Sigma operator for SBL/CSBL
%
%   Sigma = makeSigmaOp(A, alpha, beta)
%
%   INPUT:
%     A      - (M x N) sensing matrix
%     alpha  - (N x 1) ARD hyperparameters
%     beta   - noise precision
%
%   OUTPUT:
%     Sigma  - structure with operator handles:
%              .apply(v)   -> Sigma * v
%              .diag()     -> diag(Sigma)
%              .cols(J)    -> Sigma(:,J)
%              .toDense(b) -> block-wise dense reconstruction
%
%   NOTE:
%     Sigma = (beta*A'*A + diag(alpha))^(-1)
%     Using Woodbury identity:
%       Sigma = D^(-1) - D^(-1)A'(I/beta + A D^(-1) A')^(-1) A D^(-1)
%     with D = diag(alpha).

    [M, ~] = size(A);
    d  = 1 ./ alpha;              % diagonal of D^(-1)
    Bd = A .* d.';                % A*diag(d), size M x N

    % K = I/beta + A D^(-1) A' = I/beta + Bd * A'
    K  = (1/beta)*eye(M) + Bd * A.';
    Lk = chol(K, 'lower');        % Cholesky factor

    Sigma.apply   = @(v) applySigma(v, A, d, Lk);
    Sigma.diag    = @()   diagSigma(A, d, Lk, Bd);
    Sigma.cols    = @(J)  SigmaCols(J, A, d, Lk);
    Sigma.toDense = @(blk) SigmaDense(A, d, Lk, blk);
end

function y = applySigma(v, A, d, Lk)
% applySigma  Multiply Sigma by vector(s): y = Sigma * v
    W  = bsxfun(@times, d, v);        % D^(-1) * v
    AW = A * W;
    Z  = Lk \ AW; Z = Lk' \ Z;        % K^(-1) * (A*W)
    ATZ= A.' * Z;
    y  = W - bsxfun(@times, d, ATZ);
end

function ds = diagSigma(~, d, Lk, Bd)
% diagSigma  Compute diag(Sigma) efficiently
    T  = Lk \ Bd; T = Lk' \ T;        % K^(-1) * Bd
    ds = d - sum(Bd .* T, 1).';       % N x 1
end

function Scols = SigmaCols(J, A, d, Lk)
% SigmaCols  Extract selected columns Sigma(:,J)
    EJ  = zeros(length(d), numel(J));
    EJ(sub2ind(size(EJ), J(:), (1:numel(J)).')) = 1;
    WJ  = bsxfun(@times, d, EJ);
    AWJ = A * WJ;
    ZJ  = Lk \ AWJ; ZJ = Lk' \ ZJ;
    ATZJ= A.' * ZJ;
    Scols = WJ - bsxfun(@times, d, ATZJ);
end

function S = SigmaDense(A, d, Lk, blk)
% SigmaDense  Reconstruct dense Sigma (block-wise, memory expensive)
%   blk - block size (number of columns per iteration)
    N = length(d);
    if nargin<4 || isempty(blk), blk = 1024; end
    S = zeros(N, N, 'like', d);
    for j = 1:blk:N
        J = j : min(N, j+blk-1);
        S(:, J) = SigmaCols(J, A, d, Lk);
    end
end
