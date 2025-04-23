%% Choose method to generate sparse dictionary
function [phi, phi_rt] = generate_phi(method, map)
    switch(method)
        case "idw"
            [phi, phi_rt] = phi_idw(map);
        otherwise
            error('Unsupported method in generating φ: %s', method);
    end
end

% IDW interpolation
function [phi, phi_rt] = phi_idw(map)
    % Basic param
    Tx = map.Tx;
    M = width(Tx);
    N = map.size;

    % Initialize phi_rt matrix
    phi_rt = sparse(N, N);
    for i = 1:M
        phi_rt(:, Tx(i).pos) = Tx(i).gain_ln(:);
    end
    phi = phi_rt;

    % Perform IDW interpolation
    p = 2; % Set IDW parameter
    for i = 1:N
        [~, x_pos, values] = find(phi_rt(i, :));
        y_pos = i * ones(M, 1);
%         x_zero = setdiff(1:N, x_pos);
%         y_zero = i * ones(size(x_zero));
        x_zero = x_pos;
        y_zero = y_pos;

        % Perform IDW interpolation
        coordinates_zero = [x_zero(:), y_zero(:)];
        phi_temp = idw_interpolation(x_pos, y_pos, values, coordinates_zero(:, 1), coordinates_zero(:, 2), p);
        phi = phi + sparse(i, x_zero, phi_temp, N, N);
    end
end

