% Randomly select sample locations
function psi = psi_random(M, N)
% psi_random  随机从 N 个位置里不重复地选 M 个位置，每行恰好一个 1
%   M: 行数（观测数，M ≤ N）
%   N: 列数（总位置数）
    if M > N
        error('M must be <= N when sampling without replacement');
    end

    rows = (1:M)';            % 每个 1 对应的行索引
    cols = randperm(N, M)';   % 从 [1..N] 中不重复抽取 M 个列索引
    vals = ones(M,1);         % 全部赋 1

    % 一次性生成稀疏矩阵
    psi = sparse(rows, cols, vals, M, N);
    
end

