%% Choose method to generate measurement matrix
function [psi] = generate_psi(method, map, sample_rate)
    M = width(map.Tx);
    N = map.size;
    switch(method)
        case "random"
            psi = psi_random(M, N, sample_rate);
        case "mmi"
            psi = psi_mmi(M, N, sample_rate);
        otherwise
            error('Unsupported method in generating ψ: %s', method);
    end
end

% Randomly select sample locations
function psi = psi_random(M, N, sample_rate)
    % M: psi height
    % N: psi width
    psi = zeros(M, N);
    
    samples = round(sample_rate * N);
    if samples > M * N
        error('target_ones cannot be greater than the total number of elements in the matrix.');
    end
    
    for i = 1:M
        ones_positions = randperm(N, samples);  % Randomly select position
        psi(i, ones_positions) = 1;
    end
end

% Maximize mutual information
function psi = psi_mmi(M, N, sample_rate)
    % M: psi height
    % N: psi width
    psi = zeros(M, N);

    error('method psi_mmi is not available.');

end