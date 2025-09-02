function psi = psi_mmi(M, N, phi)
% Maximize mutual information (efficient & stable)
% phi: N x K  (N=候选位置数，K=字典列维)
% 返回: psi (M x N) 选点矩阵（每行 one-hot）
    % --- 兼容性：若传入的 N 与 phi 行数不一致，以 phi 行数为准 ---
    if size(phi,1) ~= N
        N = size(phi,1);
    end

    % 建议用稀疏矩阵，后续 psi*phi 依然可用且更省内存
    psi = sparse(M, N);
    selected = false(N,1);

    % 预计算：每一行的二范数平方 ||phi_i||^2
    rowNorm2 = sum(phi.^2, 2);        % N x 1

    % 已选子空间的正交基 Q（K x t）
    Q = zeros(size(phi,2), 0);

    % 所有行在 span(Q) 上的投影范数平方：proj2(i) = ||phi_i Q||^2
    proj2 = zeros(N,1);

    tol = 1e-12;   % 正交化阈值，防止病态/数值抖动
    tStart = tic;
    fprintf('MMI running ...\n');

    for t = 1:min(M, N)
        % 只在未选候选上计算增益：g_i = ||phi_i||^2 - ||phi_i Q||^2
        candidates = find(~selected);
        gains = rowNorm2(candidates) - proj2(candidates);
        gains(gains < 0) = 0;         % 数值安全：避免-1e-16这类负小数

        % 选择增益最大的候选
        [~, idx_rel] = max(gains);
        i_best = candidates(idx_rel);

        % 记录到 psi（one-hot）
        selected(i_best) = true;
        psi(t, i_best) = 1;

        % === 增量更新正交基 Q 与 proj2 ===
        % 取被选行向量并对 Q 做正交化：r_perp = r - Q*(Q'*r)
        r = phi(i_best,:).';                 % K x 1
        if ~isempty(Q)
            r = r - Q*(Q.'*r);
        end
        nr = norm(r);

        % 若有新的独立方向，则扩展 Q，并增量更新所有行的投影范数
        if nr > tol
            q_new = r / nr;                  % K x 1（单位向量）
            Q = [Q, q_new];                  % 扩展正交基
            v = phi * q_new;                 % N x 1，每行在新方向上的投影
            proj2 = proj2 + v.^2;            % 累加：||phi_i Q||^2 += (phi_i·q_new)^2
        end
    end

    fprintf('MMI duration: %.4f s\n', toc(tStart));
end
