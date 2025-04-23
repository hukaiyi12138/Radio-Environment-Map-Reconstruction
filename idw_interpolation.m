%% IDW interpolation
function interpolated_values = idw_interpolation(x, y, values, xi, yi, p)
    if nargin < 6
        p = 2;  % Default power
    end

    % Ensure x and y are column vectors and have the same length
    if length(x) ~= length(y)
        error('x and y must have the same length.');
    end
    
    % Make sure x and y are column vectors
    x = x(:);
    y = y(:);

    % 构建KDTree用于快速查找最近邻
    tree = KDTreeSearcher([x, y]);

    % 查询每个插值点的最近邻
    [idx, distances] = knnsearch(tree, [xi(:), yi(:)], 'K', length(values));

    % 防止距离为零（即插值点与已知点重合）导致除零错误
    distances(distances == 0) = 1e-10;

    % 计算权重
    weights = 1 ./ distances.^p;
    weights = weights ./ sum(weights, 2); % 将沿着行方向进行求和

    % 进行插值
    interpolated_values = sum(weights .* values(idx), 2);
end

