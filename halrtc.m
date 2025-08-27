function [X_out] = halrtc(X_in, mask)
    %—— 稀疏输入处理 ——%
    data = full(X_in);
    fprintf("rank = %d\n", rank(data)); % input rank

    %—— 二维矩阵处理 ——%
    [I, J] = size(data);
    data = reshape(data, [I, J, 1]); 
    
    %—— 掩码生成 ——%
    valid_mask = (data ~= 0);  % 找到不需要插值的位置
    valid_mask(mask) = 1; % 找到特定的不需要插值的位置

    S = double(valid_mask);  
    X = S .* data;    

    %—— 参数设置 ——%
    alpha = [1/2, 1/2];
    rho = 1e-2;
    MaxIter = 3000;
    tol = 1e-6;

    %—— 预分配内存 ——%
    X_hat = X;
    Y_cell = {zeros(I,J), zeros(I,J)};
    convergence = zeros(MaxIter,1);

    %—— 核心迭代 ——%
    % 在迭代前添加低秩增强项
    reg = 1e-5 * eye(size(X_hat(:,:,1)));  % 正则化矩阵
    tic;
    fprintf("HaLRTC running...\n");
    for it = 1:MaxIter
        fprintf("Iter = %d\n", it);
        B_cell = cell(1,2); 

        % 并行处理行/列方向
        parfor idx = 1:2  
            X_unfold = data_unfold(X_hat, idx);
            Y_unfold = Y_cell{idx};
            M = gpuArray(X_unfold + Y_unfold/rho);
            B_cell{idx} = data_refold(SVT(M, alpha(idx)/rho), idx, I, J);
        end

        % 合并结果
        avgB = (B_cell{1} + B_cell{2} - (Y_cell{1}+Y_cell{2})/rho) / 2 + reg;
        X_hat = avgB; 
        X_hat(valid_mask) = X(valid_mask);

        % 更新乘子
        Y_cell{1} = Y_cell{1} - rho*(B_cell{1} - X_hat);
        Y_cell{2} = Y_cell{2} - rho*(B_cell{2} - X_hat);

        % 更新rho与tol
        if it > 20 && mod(it,10)==0
            rho = rho * 0.935;  % 逐步衰减rho，避免梯度爆炸
            % tol = max(tol*1.01, 1e-5);
        end

        % 收敛检测
        convergence(it) = norm(X_hat(:) - avgB(:));
        if it > 10 && abs(convergence(it)-convergence(it-1)) < tol
            break; 
        end
    end

    %—— 稀疏输出 ——%
    X_out = sparse(X_hat(:,:,1));
    fprintf('HaLRTC completed in %.2f sec (%d iters)\n', toc, it);
end

%%—— 展开函数 ——%
function mat = data_unfold(tensor, mode)
    [I, J, ~] = size(tensor);
    switch mode
        case 1  % 行方向展开：I x J 矩阵
            mat = tensor(:,:,1);
        case 2  % 列方向展开：J x I 矩阵（转置）
            mat = tensor(:,:,1)';
        otherwise
            error('仅支持 mode=1 或 2');
    end
end

%%—— 折叠函数 ——%
function tensor = data_refold(mat, mode, I, J)
    switch mode
        case 1  % 行方向还原
            tensor = reshape(mat, [I, J, 1]);
        case 2  % 列方向还原（需转置回原始方向）
            tensor = reshape(mat', [I, J, 1]);
        otherwise
            error('仅支持 mode=1 或 2');
    end
end

%%—— 奇异值阈值函数——%
function Xd = SVT(X, tau)
    X_gpu = gpuArray(X);
    [U, S, V] = svd(X_gpu, 'econ');
    s = diag(S);
    s(s < 1e-10) = 1e-10;  % 防止奇异值归零
    thr = max(s - tau, 0);
    Xd_gpu = U * diag(thr) * V';
    Xd = gather(Xd_gpu);
end