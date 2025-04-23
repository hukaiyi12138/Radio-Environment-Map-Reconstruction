%% Recover signal through different methods
function [omega_est, others] = recover_signal(method, y, Phi, sigma)
    others = [];
    iter_omp = 30;
    iter_sbl = 50;
    sparsity = height(Phi);
    switch(method)
        case "omp"
            % OMP
            [omega_est, ~, ~] = omp(y, Phi, iter_omp);

        case "sbl"
            % SBL
            [omega_est, ~, ~] = sbl(y, Phi, iter_sbl, sigma);

        case "csbl"
            % CSBL & Truncation
            [omega_est_sbl, ~, ~] = csbl(y, Phi, iter_sbl, sigma);
            omega_est_trun = thresholding_sbl(omega_est_sbl); % Truncation
            omega_est = omega_est_trun; % Result

        case "msbl"
            % SBL & Truncation & MMD
            [omega_est_sbl, ~, ~] = sbl(y, Phi, iter_sbl, sigma);
            omega_est_trun = thresholding_sbl(omega_est_sbl); % Truncation
            omega_est_mmd = mmd_cluster(omega_est_trun, sparsity); % MMD method
            omega_est = omega_est_mmd; % Result

        case "cmsbl"
            % CSBL & Truncation & MMD
            [omega_est_sbl, ~, ~] = csbl(y, Phi, iter_sbl, sigma);
            omega_est_trun = truncate(omega_est_sbl); % Truncation
            omega_est_mmd = mmd_cluster(omega_est_trun, sparsity); % MMD method
            omega_est = omega_est_mmd; % Result

        otherwise
            error('Unsupported method in recover signal: %s', method);
    end

    % Trans to sparse matrix
    omega_est = sparse(omega_est);

end

% Perform truncation
function omega_est_trun = truncate(omega_est)
    % Sparse threshold
    sparse_thre = -20;
    
    % Truncation
    omega_est_trun = omega_est;
    max_omega = max(omega_est_trun);
    for i = 1:length(omega_est_trun)
        condition = 20 * log10(omega_est_trun(i) / max_omega);  % Condition
        if condition < sparse_thre 
            omega_est_trun(i) = 0;
        end
    end

end

% Perform MMD
function omega_est_mmd = mmd_cluster(omega_est, K)
        [centroids, cluster_idx] = mmd(omega_est, K);
        omega_est_mmd = centroids(cluster_idx);
end

