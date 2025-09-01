function [X_out] = halrtc(X_in, mask)
    % ===== 输入与观测集合 =====
    data = double(full(X_in));
    [I, J] = size(data);
    Omega = (data ~= 0);                 
    Omega(mask) = 1;                     % 线性索引或同形状 logical 都可

    % ====== 归一化：按“观测量级”缩放到 O(1) ======
    X_hat = zeros(I,J); 
    if any(Omega(:))
        scale = max(abs(data(Omega)));
    else
        scale = max(abs(data(:)));
    end
    data_s = data / scale;               % 归一化后的观测数据
    X_hat(Omega) = data_s(Omega);        % 初值满足约束（用归一化值）

    % ===== 参数 =====
    alpha = [0.5, 0.5];
    rho   = 1e-2;                        
    MaxIter = 1000;
    tol = 1e-5;

    % ===== 对偶变量：与“展开后的形状”一致 =====
    Y1 = zeros(I, J);                     % mode-1：I x J
    Y2 = zeros(J, I);                     % mode-2：J x I

    % ===== 迭代 =====
    fprintf("HaLRTC running...\n");
    t0 = tic;
    M1_prev = zeros(I,J); M2_prev = zeros(J,I);

    for it = 1:MaxIter
        % --- 展开（在归一化尺度上）---
        X1 = data_unfold_2d(X_hat, 1);    % I x J
        X2 = data_unfold_2d(X_hat, 2);    % J x I

        % --- M 步（每个模做 SVT，阈值 = alpha_i / rho）---
        M1 = SVT_cpu(X1 + Y1/rho, alpha(1)/rho);
        M2 = SVT_cpu(X2 + Y2/rho, alpha(2)/rho);

        % --- X 步（合并 + 投影到观测约束）---
        X_bar = ( ...
            data_refold_2d(M1 - Y1/rho, 1, I, J) + ...
            data_refold_2d(M2 - Y2/rho, 2, I, J) ) / 2;

        X_hat = X_bar;
        X_hat(Omega) = data_s(Omega);     % 硬投影
        X_hat(~Omega) = min(max(X_hat(~Omega), 0), 1);

        % --- 对偶更新：在与展开一致的空间里 ---
        X1_new = data_unfold_2d(X_hat, 1);        % I x J
        X2_new = data_unfold_2d(X_hat, 2);        % J x I
        Y1 = Y1 + rho * (X1_new - M1);
        Y2 = Y2 + rho * (X2_new - M2);

        % --- 收敛判据（原始/对偶残差）---
        r1 = norm(X1_new - M1, 'fro'); 
        r2 = norm(X2_new - M2, 'fro');
        r  = sqrt(r1^2 + r2^2);                   % primal residual

        s1 = rho * norm(M1 - M1_prev, 'fro');
        s2 = rho * norm(M2 - M2_prev, 'fro');
        s  = sqrt(s1^2 + s2^2);                   % dual residual

        M1_prev = M1; M2_prev = M2;

        if r < tol && s < tol
            fprintf('Converged at iter %d | r=%.3e, s=%.3e\n', it, r, s);
            break;
        end
        if mod(it,50)==0
            fprintf('iter %d | r=%.3e, s=%.3e\n', it, r, s);
        end
    end

    % ===== 反缩放输出 =====
    X_out = sparse(X_hat * scale);
    fprintf('HaLRTC completed in %.2fs (%d iters)\n', toc(t0), it);
end

% ===== 仅矩阵情形下的 unfold/refold =====
function mat = data_unfold_2d(X, mode)
    switch mode
        case 1, mat = X;      
        case 2, mat = X.';    
        otherwise, error('mode must be 1 or 2');
    end
end

function X = data_refold_2d(mat, mode, I, J)
    switch mode
        case 1, X = reshape(mat, [I, J]);
        case 2, X = reshape(mat.', [I, J]);
        otherwise, error('mode must be 1 or 2');
    end
end

% ===== 稳定的 CPU 版 SVT =====
function Xd = SVT_cpu(X, tau)
    [U,S,V] = svd(X, 'econ');
    s = diag(S);
    s = max(s - tau, 0);
    if isempty(s)
        Xd = zeros(size(X));
    else
        Xd = U * diag(s) * V';
    end
end
