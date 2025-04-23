%% Maximize mutual information algorithm
function [centroids_out, cluster_idx] = mmd(data_in, K)
    % K: 簇的数量
    % centroids: 簇中心
    % cluster_idx: 每个数据点所属的簇的索引

    % 数据大小
    n = length(data_in);

    % Step 0: 预处理，标准化
    max_data = max(data_in);
    min_data = min(data_in);
    data = (data_in - min_data) / (max_data - min_data); % 标准化至区间[0,1]

    % Step 1: 初始化K个簇中心，随机选择K个数据点作为初始簇中心
    random_idx = randperm(n, K); % 随机选择K个数据点
    centroids = data(random_idx); % 初始簇中心为随机选择的点

    % Step 2: 初始化隶属度矩阵
    dist = zeros(n, K); % 距离矩阵，存储每个数据点与每个簇中心的距离
    cluster_idx = zeros(n, 1); % 每个数据点所属的簇索引

    max_iter = 10; % 最大迭代次数
    tol = 1e-6; % 收敛阈值
    prev_centroids = centroids; % 上一次的簇中心

    for iter = 1:max_iter
        % Step 3: 计算每个数据点与每个簇中心的距离
        for i = 1:n
            for j = 1:K
                % 计算数据点与簇中心之间的距离，避免除以零
                dist(i, j) = abs(data(i) - centroids(j));
            end
        end
        
        % Step 4: 根据最小距离分配每个数据点到最接近的簇
        [~, cluster_idx] = min(dist, [], 2); % 找到最小距离的簇中心
        
        % Step 5: 更新簇中心
        for j = 1:K
            % 获取当前簇内的所有数据点
            points_in_cluster = data(cluster_idx == j);
            
            % 如果该簇没有数据点，跳过更新或给簇中心一个默认值
            if isempty(points_in_cluster)
                % 可以选择保留原簇中心不变，或者选择重新随机选取一个点
                continue;
            else
                % 计算簇内数据点的均值作为新的簇中心
                centroids(j) = mean(points_in_cluster);
            end
        end
        
        % Step 6: 检查收敛条件
        if max(abs(centroids - prev_centroids)) < tol
            break;
        end
        prev_centroids = centroids; % 更新前一轮的簇中心
    end

    % 还原至原始范围
    centroids_out = centroids * (max_data - min_data) + min_data;

end
