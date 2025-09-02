function phi = kriging(phi_rt, map)
%KRIGING  Ordinary Kriging for each selected Tx column of phi_rt.
%   phi = KRIGING(phi_rt, map)
%   - phi_rt: N x N 稀疏矩阵，已在 partPos 有值、在 interPos 为 0（见 generate_phirt）
%   - map:    需要用到 Nx, Ny, size, partPos, interPos, build, selectedTxPos
%   - 返回:   N x N（稀疏）矩阵，仅 selectedTxPos 列被填充

    N   = map.size;
    Nx  = map.Nx;
    %Ny = map.Ny;  % 如需用到可打开
    cols = map.selectedTxPos(:);

    % 参数（可按需要微调）
    nNeighbors    = 50;                      % 每个目标点使用的近邻个数（20~80 之间比较常用）
    rangeFactor   = 0.25;                    % 相关长度 = rangeFactor * max(Nx, Ny)
    nuggetFactor  = 0.05;                    % nugget = nuggetFactor * var(z)
    jitter        = 1e-10;                   % 数值稳定性

    % 预计算所有索引的 (x,y)
    [xAll, yAll] = idx2xy( (1:N).', Nx );

    obsIdx  = map.partPos(:);                % 已观测
    targIdx = map.interPos(:);               % 待估（不包含建筑）
    % 若你希望跳过建筑格子，确保 interPos 本身就不含 build

    phi = phi_rt;                            % 在此基础上填充
    % 对每个被选中的发射机列做克里金
    for cc = 1:numel(cols)
        j = cols(cc);

        % 观测 z，注意把稀疏取为 full
        z_obs = full(phi_rt(obsIdx, j));
        % 若没有有效观测（极端情况），跳过
        if nnz(z_obs) < 3
            warning('Kriging: column %d has too few observations, skip.', j);
            continue;
        end

        % 观测点坐标
        Xobs = [xAll(obsIdx), yAll(obsIdx)];

        % 经验尺度参数
        sig2   = var(z_obs);                         % partial sill 估计
        if sig2 <= eps
            % 几乎常数场：全图直接用常数
            phi(targIdx, j) = z_obs(1);
            continue;
        end
        range  = rangeFactor * max(map.Nx, map.Ny);  % 相关长度
        nugget = nuggetFactor * sig2;

        % 逐个目标点做局部 OK，使用 KNN 近邻
        for t = 1:numel(targIdx)
            id0 = targIdx(t);
            x0  = xAll(id0); y0 = yAll(id0);

            % 与所有观测点的距离（欧氏）
            d = hypot(Xobs(:,1) - x0, Xobs(:,2) - y0);

            % 选最近的 nNeighbors 个观测
            k = min(nNeighbors, numel(d));
            [dSort, nnIdx] = mink(d, k);
            Xnn = Xobs(nnIdx, :);
            znn = z_obs(nnIdx);

            % 如果邻居方差过小或邻居过少，退化为 IDW 兜底
            if var(znn) <= eps || k < 3
                w = 1 ./ max(dSort, 1e-6).^2;
                phi(id0, j) = sum(w .* znn) / sum(w);
                continue;
            end

            % --- Ordinary Kriging（协方差形式）---
            % 协方差模型：C(h) = sig2 * exp(-h/range)
            Hnn = pairwise_dist(Xnn, Xnn);
            Sigma = cov_exp(Hnn, sig2, range);

            % nugget 加在对角上；再加一点 jitter 防止奇异
            Sigma(1:k+1:end) = Sigma(1:k+1:end) + (nugget + jitter);

            % 目标与邻居的协方差向量
            h0 = hypot(Xnn(:,1) - x0, Xnn(:,2) - y0);
            c0 = cov_exp(h0, sig2, range);

            % OK 线性系统：
            % [Sigma  1] [w]   = [c0]
            % [ 1^T   0] [mu]    [ 1]
            A = [Sigma,            ones(k,1); 
                 ones(1, k),       0        ];
            b = [c0; 1];

            % 解权重
            wl = A \ b;
            w  = wl(1:k);

            % 预测
            zhat = w' * znn;

            % 物理上 phi ≥ 0（线性功率），做个下界截断
            if zhat < 0, zhat = 0; end
            phi(id0, j) = zhat;
        end
    end

    % 建筑位置保持 0（如需）
    if isfield(map, 'build') && ~isempty(map.build)
        % phi(map.build, :) = 0;  % 通常 interPos/partPos 已经排除建筑，无需重复处理
    end

    % 稀疏化（可选）
    phi = sparse(phi);
end

% ------ helpers ------

function C = cov_exp(h, sig2, range)
% 指数核协方差：C(h) = sig2 * exp(-h / range)
    C = sig2 .* exp(-h ./ max(range, eps));
end

function D = pairwise_dist(X, Y)
% X: m×2, Y: n×2 -> D: m×n 欧氏距离
    Dx = X(:,1) - Y(:,1).';
    Dy = X(:,2) - Y(:,2).';
    D  = hypot(Dx, Dy);
end

function [x, y] = idx2xy(idx, Nx)
% 把线性索引（1..N）转为 (x,y)（与 itc 等价）
    x = mod(idx-1, Nx) + 1;
    y = floor((idx-1) / Nx) + 1;
end