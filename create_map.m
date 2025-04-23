%% Generate source
direct_map = "source";
file_map = "map.mat";
if ~exist(direct_map, "dir")
    mkdir(direct_map);
end
file_map_path = sprintf("%s/%s", direct_map, file_map);

sparse_dict = struct("sparsity", [],"method", [], "phi", [], "phi_rt", []);

if ~exist(file_map_path,"file")     
    % Generate map within different sparsity
    K = [4, 8, 12, 16];
    for i = 1:length(K)
        map(i) = generate_map_2D(K(i));

        % Generate phi according to map
        sparse_dict(i).sparsity = K(i);
        [sparse_dict(i).phi, sparse_dict(i).phi_rt] = generate_phi(method_phi, map(i));

    end

    % Generate measurement matrix
    sample_rate = linspace(0.01, 0.15, 15);
    for j = 1:length(sample_rate)
        [psi] = generate_psi(method, map, sample_rate);
    end


    save(file_map, "map");

end



