%% Use real map dataset
clear;

% Set params
% K = 4;                     % {4, 8, 12, 16}
sample_rate = 0.05;         % 0.01 ~ 0.15
sigma = 0.1;                % Noise power: sigma = σ^2
method_phi = "idw";         % {idw}
method_psi = "random";      % {random, mmi}
method_recov = "sbl";     % {omp, sbl, csbl, msbl, cmsbl}
% Nx = 50;
% Ny = 50;

% Load map file
direct_input = "dataset";
map_name = "";


% [map] = generate_map2D(K, Nx, Ny);
omega_real = map.omega_real;

% Generate sparse dictionary
[phi, phi_rt] = generate_phi(method_phi, map);

% Generate measurement matrix
[psi] = psi_random(numel(map.Tx), map.size, sample_rate);

% Transmit process
Phi = psi * phi; % Sensing matrix
y = Phi * omega_real; % Observation vector

% Recover signal
[omega_est] = recover_signal(method_recov, y, Phi, sigma);

% Evaluation
mse = norm(omega_real - omega_est) / norm(omega_real);
mse_db = 10 * log10(mse);

% Save result
direct_output = "result";
if ~exist(direct_output,"dir")
    mkdir(direct_output);
end
result_name = sprintf('%s/K=%d_r=%.2f_%s_%s.mat', direct_output, K, sample_rate, method_psi, method_recov);
save(result_name);