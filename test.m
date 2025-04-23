%% test: run a simple experiment
% Set trial params
K = 4;
sample_rate = 0.05;
sigma = 0.1; % Noise power: sigma = σ^2
method_phi = "idw";
method_psi = "random"; % {random, mmi}
method_recov = "cmsbl"; % {omp, sbl, csbl, msbl, cmsbl}
Nx = 50;
Ny = 50;

% Generate map and Tx info
[map] = generate_map2D(K, Nx, Ny);
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
direct_name = "test";
if ~exist(direct_name,"dir")
    mkdir(direct_name);
end
result_name = sprintf('%s/K=%d_r=%.2f_%s_%s.mat',direct_name, K, sample_rate, method_psi, method_recov);
save(result_name);
