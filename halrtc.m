function [X_hat_sparse] = halrtc(X_in)
    %—— 把任何稀疏输入都转成密集 ——%
    if issparse(X_in)
        data = full(X_in);
    else
        data = X_in;
    end

    %—— 参数 ——%
    alpha1 = 1/3; alpha2 = 1/3; alpha3 = 1/3;
    rho     = 1e-2;
    MaxIter = 100;
    
    %—— 初始化 ——%
%     [I,J,K] = size(data);
    d = ndims(data);
    if d == 2
        [I,J] = size(data);
        K = 1;
        data = reshape(data, [I, J, 1]);
    elseif d == 3
        [I,J,K] = size(data);
    else
        error('halrtc: only supports 2D or 3D input');
    end
    
    allZeroCols = all(data==0, 1);
    maskUnobs = (data==0) & ~allZeroCols;
    S = ones(size(data));
    S(maskUnobs) = 0;

    obs     = data;
%     S       = round(rand(I,J,K) + 0.3);
    X       = S .* obs;
    
    X_hat = X;

    Y1 = zeros(I,J,K);  Y2 = Y1;  Y3 = Y1;
    convergence = zeros(MaxIter,1);
    
    %—— 迭代 ——%
    tic;
    fprintf("HaLRTC interpolating ...\n");
    for it = 1:MaxIter
        % fprintf("HaLRTC iter = %d\n", it);
        B1 = t_refold( SVT( t_unfold(X_hat,1) + t_unfold(Y1,1)/rho, alpha1/rho ), 1, I,J,K );
        B2 = t_refold( SVT( t_unfold(X_hat,2) + t_unfold(Y2,2)/rho, alpha2/rho ), 2, I,J,K );
        B3 = t_refold( SVT( t_unfold(X_hat,3) + t_unfold(Y3,3)/rho, alpha3/rho ), 3, I,J,K );
        
        avgB = (B1+B2+B3 - (Y1+Y2+Y3)/rho) / 3;
        X_hat = (1-S).*avgB + S.*X;

        Y1 = Y1 - rho*(B1 - X_hat);
        Y2 = Y2 - rho*(B2 - X_hat);
        Y3 = Y3 - rho*(B3 - X_hat);
        
        convergence(it) = sum( X_hat(maskUnobs).^2 );
    end
    
    %—— 输出：如果是 3D，就做 cell array；否则直接 sparse 矩阵 ——%
    if K > 1
        X_hat_sparse = cell(1,K);
        for k = 1:K
            X_hat_sparse{k} = sparse( X_hat(:,:,k) );
        end
    else
        X_hat_sparse = sparse( X_hat(:,:,1) );
    end
    
    % 计算 RMSE 和 MRE
    err = obs(maskUnobs) - X_hat(maskUnobs);
    RMSE = sqrt(mean(err.^2));
    MRE  = sum(abs(err)) / sum(X_hat(maskUnobs));
    
    fprintf('HaLRTC finished\n');
    fprintf("HaLRTC duration: %.4f sec/%d iters\n", toc, it);
    fprintf('RMSE = %.4f  (veh/15min)\n', RMSE);
    fprintf('MRE  = %.4f\n', MRE);

end

%%——— 函数：张量展开 ———%%
function mat = t_unfold(tensor, mode)
    [I, J, K] = size(tensor);
    switch mode
        case 1
            % 沿 mode-1：得到 I × (J*K)
            mat = reshape(permute(tensor, [2 3 1]), J, K*I).';
        case 2
            % 沿 mode-2：得到 J × (I*K)
            mat = reshape(permute(tensor, [1 3 2]), I, K*J).';
        case 3
            % 沿 mode-3：得到 K × (I*J)
            mat = reshape(permute(tensor, [1 2 3]), I, J*K).';
        otherwise
            error('t_unfold: mode must be 1,2,3');
    end
end

%%——— 函数：张量折叠 ———%%
function tensor = t_refold(mat, mode, I, J, K)
    switch mode
        case 1
            tmp = reshape(mat.', [J, K, I]);
            tensor = permute(tmp, [3 1 2]);
        case 2
            tmp = reshape(mat.', [I, K, J]);
            tensor = permute(tmp, [1 3 2]);
        case 3
            tmp = reshape(mat.', [I, J, K]);
            tensor = tmp;
        otherwise
            error('t_refold: mode must be 1,2,3');
    end
end

%%——— 函数：奇异值阈值化 SVT ———%%
function Xd = SVT(X, tau)
    [U, S, V] = svd(X, 'econ');
    S_threshold = max(diag(S) - tau, 0);
    Xd = U * diag(S_threshold) * V';
end

