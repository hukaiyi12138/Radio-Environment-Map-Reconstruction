%% Choose method to generate measurement matrix
function [psi] = generate_psi(method, map, sample_rate, phi)
    M = width(map.Tx);
    N = map.size;
    switch(method)
        case "random"
            psi = psi_random(M, N, sample_rate);
        case "mmi"
            psi = psi_mmi(M, N, sample_rate, phi);
        otherwise
            error('Unsupported method in generating ψ: %s', method);
    end
end
