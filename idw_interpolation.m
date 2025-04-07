%% IDW
% % 假设 image_samples 是一个二维矩阵
% % 提取非零值的位置和对应的值
% [x, y] = find(image_samples);  % 找到非零值的位置
% values = image_samples(sub2ind(size(image_samples), x, y));  % 对应的值
% 
% % 找出0值的位置
% % 所有位置的坐标
% [all_x, all_y] = meshgrid(1:size(image_samples, 1), 1:size(image_samples, 2));
% all_coords = [all_x(:), all_y(:)];
% 
% % 已采样点的坐标
% sampled_coords = [x, y];
% 
% % 找到剩余的未采样点（0值位置）
% remaining_coords = setdiff(all_coords, sampled_coords, 'rows');
% zero_x = remaining_coords(:, 1);  % 未采样点的 x 坐标
% zero_y = remaining_coords(:, 2);  % 未采样点的 y 坐标
% 
% % 使用IDW进行插值
% interpolated_values = idw_interpolation(x, y, values, zero_x, zero_y);

function interpolated_values = idw_interpolation(x, y, values, xi, yi, power)
    if nargin < 6
        power = 2;  % Default power
    end

    % 构建KDTree用于快速查找最近邻
    tree = KDTreeSearcher([x(:), y(:)]);

    % 查询每个插值点的最近邻
    [idx, distances] = knnsearch(tree, [xi(:), yi(:)], 'K', length(values));

    % 防止距离为零（即插值点与已知点重合）导致除零错误
    distances(distances == 0) = 1e-10;

    % 计算权重
    weights = 1 ./ distances.^power;
    weights = weights ./ sum(weights, 2); % 将沿着行方向进行求和

    % 进行插值
    interpolated_values = sum(weights .* values(idx), 2);
end
