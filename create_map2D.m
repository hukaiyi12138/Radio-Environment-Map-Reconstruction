%% Generate map source
clear;
Nx = 50;
Ny = 50;
Klist = [4, 8, 12, 16]; % sparsity set

% set path
direct_map = "test_map";
file_map = "map.mat";
if ~exist(direct_map,"dir")
    mkdir(direct_map);
end
file_map_path = fullfile(direct_map, file_map);

% load all the maps
if isfile(file_map_path)
    S = load(file_map_path, 'map');
    map = S.map;    % struct array with fields: K, Nx, Ny, grid
else
    map = struct('sparsity',{},'Nx',{},'Ny',{},'grid',{});
end

% check if map exists
didAppend = false;
for k = Klist
    isExist = any( [map.sparsity]==k & [map.Nx]==Nx & [map.Ny]==Ny );
    if ~isExist
        fprintf("  Generating map for K=%d, Nx=%d, Ny=%d ...\n", k, Nx, Ny);
        G = generate_map2D(k, Nx, Ny);
        newEntry = struct('sparsity', k, 'Nx', Nx, 'Ny', Ny, 'grid', G);
        map(end+1) = newEntry;
        didAppend = true;
    else
        fprintf("Skip map for K=%d, Nx=%d, Ny=%d (already exists)\n", k, Nx, Ny);
    end
end

% save new map data
if didAppend
    save(file_map_path, 'map', '-v7.3');
    fprintf("Saved updated mapDB (%d entries total)\n", numel(map));
else
    fprintf("No new entries — mapDB untouched.\n");
end

% close
clear;
