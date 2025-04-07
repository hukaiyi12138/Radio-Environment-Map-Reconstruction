%% Generate map information
function generate_map(map_index, x, y, dx, dy)
    fprintf('map = %d, size = %d*%d, start_point=(%d,%d)\n', map_index, dx, dy, x, y);

    % Ensure the 'output' folder exists
    if ~exist('output', 'dir')
        mkdir('output');
    end

    % Create a folder for saving images if it does not already exist
    output_folder = sprintf('output/map_%d', map_index);
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Check if the .mat file exists
    mat_file = sprintf('%s/map%d.mat', output_folder, map_index);
    if exist(mat_file, 'file')
        fprintf('File map%d.mat already exists.\n', map_index);
        return;  % Exit the function if the .mat file already exists
    end

    % Initialize
    Nxx = 256;
    Nyy = 256;
    map_size = Nxx * Nyy; % Total size of the input data
    Pt = 30; % Transmitter power Pt = 30 dBm

    % Read transmitter position information
    antenna_name = sprintf('RadioMap3DSeer/antenna/%d.json', map_index);
    Tx_3Ddata = jsondecode(fileread(antenna_name)); % Transmitter position information
    Tx = table(Tx_3Ddata(:, 1), Tx_3Ddata(:, 2), (1:height(Tx_3Ddata))', 'VariableNames', {'X', 'Y', 'index'});
    Tx.Y = Nxx - Tx.Y;
    M = height(Tx); % Number of transmitters

    % Calculate gain values at each transmitter position
    gain_min = -100;
    gain_max = -30;
    gain_value = zeros(map_size, M); % Units in dB
    for k = 1:M
        gain_img = imread(sprintf('RadioMap3DSeer/gain/%d_%d.png', map_index, k-1));
        gain_value(:, k) = gain_min + (gain_max - gain_min) * (double(gain_img(:)) / 255); % Compute mapped gain values
    end

    % Calculate ideal RSS (dBm)
    Pt_mW = 10^(Pt/10);
    gain = 10 .^ (gain_value ./ 10); % Linear gain values
    Pr_total = sum(Pt_mW .* gain, 2); % Sum across rows to get total RSS at each position (mW)
    Pr_dBm = 10 * log10(Pr_total); % Convert to RSS (dBm)
    RSS_entire = reshape(Pr_dBm,Nxx,Nyy);

    % Determine ROI boundaries
    Nx1 = x;
    Nx2 = Nx1 + dx;
    Ny1 = y;
    Ny2 = Ny1 + dy;
    Nx = Nx2 - Nx1 + 1;
    Ny = Ny2 - Ny1 + 1;
    N = Nx * Ny;

    % Determine transmitters within ROI in global coordinates
    roi_cond = (Nx1 < Tx.X) & (Tx.X < Nx2) & (Ny1 < Tx.Y) & (Tx.Y < Ny2);
    roi_tx_x = Tx.X(roi_cond);
    roi_tx_y = Tx.Y(roi_cond);

    % Convert global coordinates to ROI coordinates
    Tx_roi_x = roi_tx_x - Nx1;
    Tx_roi_y = roi_tx_y - Ny1;
    Tx_roi_index = cti(Tx_roi_x, Tx_roi_y, Nx);
    Tx_roi = table(Tx_roi_x, Tx_roi_y, Tx_roi_index,'VariableNames', {'X', 'Y', 'index'});

    % Generate all one-dimensional coordinates for the ROI region using meshgrid
    [x_grid, y_grid] = meshgrid(Nx1:Nx2, Ny1:Ny2);  % Create grid coordinates
    roi = cti(x_grid(:), y_grid(:), Nxx);  % Calculate the unique one-dimensional index for each grid point
    
    % Calculate the gain values within the ROI
    gain_roi = gain(roi, Tx.index);
    Pr_roi_dBm = Pr_dBm(roi); % RSS (dBm) values within the ROI
    RSS_ideal = reshape(Pr_roi_dBm, Nx, Ny);  % Reshape to match the ROI dimensions

    % Construct sparse signal ω (N*1)
    omega_real = sparse(Tx.index, 1, Pt_mW, N, 1);

    % Initialize phi_rt matrix
    phi_rt = sparse(N, N);
    phi_rt(:, Tx.index) = gain_roi;
    phi = phi_rt;

    % Perform IDW interpolation (you may need to adjust based on your specific needs)
    p = 1; % Set IDW parameter
    tic;
    for i = 1:N
        fprintf('Generate phi: %d/%d\n', i, N);
        [~, x_pos, values] = find(phi_rt(i, :)); % Extract non-zero positions and values
        y_pos = i * ones(M, 1);
        x_zero = setdiff(1:N, x_pos);  % Remove existing non-zero positions
        y_zero = i * ones(size(x_zero));

        % Perform IDW interpolation
        coordinates_zero = [x_zero(:), y_zero(:)];
        phi_temp = idw_interpolation(x_pos, y_pos, values, coordinates_zero(:, 1), coordinates_zero(:, 2), p);
        phi = phi + sparse(i, x_zero, phi_temp, N, N);
    end
    toc;

    % Save data
    save(sprintf('%s/map%d.mat', output_folder, map_index));
    fprintf('All data has been saved to %s/map%d.mat\n', output_folder, map_index);

end

% Coordinate to index
function index = cti(x, y, Nx)
    index = (x - 1) * Nx + y;
end