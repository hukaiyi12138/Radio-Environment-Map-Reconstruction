%% Choose method to generate sparse dictionary
function [phi, phi_rt, phi0] = generate_phi(method, map)
    fprintf('\n---- Generating φ ----\n');

    % Generate phi_rt
    [phi_rt, phi0] = generate_phirt(map);

    switch(method)
        case "idw"
            % % Choose the best param p
            % switch(map.name)
            %     case "tinymap"
            %         p = 12;
            %     case "largemap"
            %         p = 11;
            %     case "hugemap"
            %         p = 10;
            %     otherwise
            %         p = 2;
            % end
            p = 2;
            [phi] = phi_idw(phi_rt, map, p);

        case "halrtc"
            % For every single tx gain
            phi = sparse(map.size, map.size);
            col = map.selectedTxPos(:);
            for i = 1:numel(col)
                mat = reshape(phi_rt(:, col(i)), map.Nx, []);
                [mat2] = halrtc(mat, map.build);
                phi(:, col(i)) = mat2(:);
            end

        case "kriging"
            [phi] = kriging(phi_rt, map);

        otherwise
            error('Unsupported method in generating φ: %s', method);
    end

    phi = sparse(phi);

end

% Generate phi_rt and phi0
function [phi_rt, phi0] = generate_phirt(map)
    N    = map.size;
    
    phi0 = sparse(N, N);
    phi0(:, map.selectedTxPos) = map.phi(:, map.selectedTxPos);

    phi_rt = phi0;
    phi_rt(map.interPos, :) = 0; % unselected position
end
