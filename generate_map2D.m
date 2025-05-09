%% Generate and save map data
function [map] = generate_map2D(K, Nx, Ny)
    % map size: Nx * Ny
    % K: number of Tx / sparsity

    % Set Tx parameters
    Pt = 2000; % Transmit power = 2W
    freq = 1e9; % Frequency = 1 GHz

    % Set map size
    N = Nx * Ny;
    map = struct('sparsity', K ,'width', Nx, 'height', Ny, 'size', N, 'Tx', [], 'omega_real', []); % Output map info

    % Set Tx information
    % gain: dB value
    % gain_ln: linear gain value
    Tx(K) = struct('power', [], 'freq', [], 'id', [], 'x', [], 'y', [], 'pos', [], 'gain', [], 'gain_ln', []);
    positions = [];

    for i = 1:K
        while true
            % Generate random coordinate
            x = randi([1, Nx]);
            y = randi([1, Ny]);
            
            % Ensure no duplicate positions
            if ~judge_pos(positions, x, y)
                positions = [positions; x, y];  % Add to the list of used positions
                break;
            end
        end
        
        % Set transmitter parameters
        Tx(i).power = Pt;
        Tx(i).freq = freq;
        Tx(i).id = i;
        Tx(i).x = x;
        Tx(i).y = y;
        Tx(i).pos = (x - 1) * Nx + y;  % Convert (x, y) to a linear index

        % Calculate path loss for every Tx in every position
        Tx(i).gain = zeros(Nx, Ny);
        Tx(i).gain_ln = zeros(Nx, Ny);

        for x = 1:Nx
            for y = 1:Ny
                % Calculate the distance from each Tx to each point (x, y)
                distance = sqrt((Tx(i).x - x)^2 + (Tx(i).y - y)^2);  % Euclidean distance
                
                % Apply the free space path loss formula (in dB)
                if distance == 0
                    PL = 0;  % No loss at the transmitter location
                else
                    PL = 20 * log10(distance) + 20 * log10(freq) - 147.55;  % Path loss in dB
                end
                Tx(i).gain(x, y) = ln_to_db(Tx(i).power) - PL;  % Received signal strength (dBm)
                Tx(i).gain_ln(x, y) = db_to_ln(Tx(i).gain(x, y)); % linear gain
            end
        end
    end
    map.Tx = Tx;

    % Generate omega_real
    omega_real = sparse([Tx.pos], ones(1, K), [Tx.power], N, 1); % size: N*1
    map.omega_real = omega_real;

end

% Judge Tx position already exists or not
function mark = judge_pos(positions, x, y)
    mark = false;
    for i = 1:size(positions, 1)
        if positions(i, 1) == x && positions(i, 2) == y 
            mark = true;
        else
            mark = false;
        end
    end
end

% Convert dB to linear
function output = db_to_ln(input)
    output = 10 .^ (input ./ 10);
end

% Convert linear to dB
function output = ln_to_db(input)
    output = 10 .* log10(input);
end

