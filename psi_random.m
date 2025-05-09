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