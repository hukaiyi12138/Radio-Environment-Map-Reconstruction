%% IDW interpolation - 2D
function phi = phi_idw(phi_rt, map, p)
    if nargin<3, p = 3; end

    N    = map.size;
    Nx   = map.Nx;
    cols = map.selectedTxPos;   % 只对这 K 列插值
    phi  = phi_rt;              % 从运行时字典开始

    % 预先算好所有格点的 (x,y)
    idx = (1:N).';
    [x_all, y_all] = itc(idx, Nx);

    fprintf("IDW interpolating with p = %d ...\n", p); 
    tic;
    for j = cols
        % 已有观测值的行 & 值
        [row_k, ~, v_k] = find(phi_rt(:, j));
        if isempty(row_k) || numel(row_k)==N
            continue
        end

        % 待插值行：只取 map.interPos 中还没观测到的
        unknown = setdiff(map.interPos, row_k);
        if isempty(unknown)
            continue
        end

        % 坐标拆分
        xk = x_all(row_k);    yk = y_all(row_k);
        xu = x_all(unknown);  yu = y_all(unknown);

        % 调用 IDW 函数
        v_u = idw_interpolation(xk, yk, v_k, xu, yu, p);

        % 写回稀疏矩阵
        phi(unknown, j) = v_u;

        fprintf("IDW: col %d done (%d pts) in %.2f s\n", j, numel(unknown), toc);
    end
    fprintf("IDW finished in %.2f s\n", toc);
end

% Convert index to coordinate - 2D
function [x, y] = itc(idx, Nx)
    x = mod(idx-1, Nx) + 1;
    y = floor((idx-1) / Nx) + 1;
end
