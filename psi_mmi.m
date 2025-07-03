% Maximize mutual information
function psi = psi_mmi(M, N, phi)
    psi = zeros(M, N);                  % 最终的测量矩阵
    selected = false(1, N);             % 已选位置掩码
    Phi_t = [];                         % 已选位置对应的字典行组成的 Φ_t
    
    tic;
    fprintf('MMI running ...\n');
    for t = 1:M
        candidates = find(~selected);  % 剩余候选位置
        if t == 1
            % 第一次选取时，只需最大化 ‖ϕ_i‖^2 
            gains = sum(phi(candidates,:).^2, 2);
        else
            % 构造投影矩阵 P = I - Φ_t'*(Φ_t*Φ_t')^{-1}*Φ_t
            A = Phi_t * Phi_t';
            P = eye(N) - Phi_t' * (A \ Phi_t);
            % 计算每个候选位置的信息增益
            gains = sum((phi(candidates,:) * P) .* phi(candidates,:), 2);
        end
        % 选取增益最大的索引
        [~, idx_rel] = max(gains);
        i_best = candidates(idx_rel);

        % 标记并更新 psi 与 Φ_t
        selected(i_best) = true;
        psi(t, i_best) = 1;
        Phi_t = [Phi_t; phi(i_best, :)];
    end
    fprintf('MMI duration: %.4lf', toc);

end
