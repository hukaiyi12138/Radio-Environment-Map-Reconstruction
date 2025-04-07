%% OMP algorithm
% y:观测向量
% A:感知矩阵
% t:迭代次数
function [theta, rn_norm, converge_point] = omp(y, A, t)
    % 设置收敛的容忍误差
    tol = 1e-3;

    % 开辟存储空间
    rn_norm = zeros(1, t);  % 存储每次迭代的r_n范数值
    converge_point = 1; % 初始化收敛点

    [y_rows, y_columns] = size(y);
    if y_rows < y_columns
        y = y';     % y应该是个列向量
    end
    [m, n] = size(A);  % A=m*n;

    theta = zeros(n, 1);  % 存储
    At = zeros(m, t);     % 用来迭代过程中存储A被选择的列
    Pos_theta = zeros(1, t);  % 用来迭代过程中存储A被选择的列序号
    r_n = y;              % 初始化残差为y

    for i = 1:t
        product = A' * r_n;  % 感知矩阵A各列与残差的内积
        [~, pos] = max(abs(product));  % 找到最大内积绝对值，即与残差最相关的列
        At(:, i) = A(:, pos);         % 存储这一列
        Pos_theta(i) = pos;           % 储存这一列的序号
        A(:, pos) = zeros(m, 1);      % 清零A的这一列，其实此行可以删除，因为其与残差正交

        % 使用lsqnonneg确保theta是非负的
        theta_ls = lsqnonneg(At(:, 1:i), y);  % 最小二乘解，确保theta是非负的
        r_n = y - At(:, 1:i) * theta_ls;  % 更新残差

        % 计算r_n的范数值
        rn_norm(i) = norm(r_n);

        % 检查收敛条件
        if i >= 2 && rn_norm(i) < rn_norm(i-1)
            converge_point = i;
        end
        if rn_norm(i) < tol
            converge_point = i;
            break;
        end
    end

    % 恢复最终的theta
    for j = 1:i  % 直到当前迭代次数为止
        theta(Pos_theta(j)) = theta_ls(j);  % 用选定位置的theta_ls值填充theta
    end

end
