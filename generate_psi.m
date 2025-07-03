%% Choose method to generate measurement matrix
function [psi] = generate_psi(method, map, sample_rate, phi)
    fprintf('\n---- Generating ψ ----\n');
    fprintf('Method: %s\n', method);

    N = map.size;
    M  = round(sample_rate * N);
    switch(method)
        case "random"
            psi = psi_random(M, N);
        case "mmi"
            psi = psi_mmi(M, N, phi);
        otherwise
            error('Unsupported method in generating ψ: %s', method);
    end
end
