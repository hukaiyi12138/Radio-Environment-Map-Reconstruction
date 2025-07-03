%% Generate 2D map
function map = generate_map2D(mapname)
    % set file directory
    switch (mapname)
        case "tinymap"
            file_tx = 'tinymap/tx.txt';
            file_rss = 'tinymap/tinymap.power.t%03d_03.r001.p2m';
            Nx = 20; Ny = 20;
            unit = 10.0;
        case "largemap"
            file_tx = 'largemap/largemap.txt';
            file_rss = 'largemap/largemap.power.t%03d_09.r005.p2m';
            Nx = 100; Ny = 100;
            unit = 100.0;
        case "hugemap"
            file_tx = 'hugemap/hugemap.txt';
            file_rss = 'hugemap/hugemap.power.t%03d_02.r003.p2m';
            Nx = 256; Ny = 256;
            unit = 256.0;
        otherwise
            error("Incorrect map name.\n");
    end

    % Target: trans to Nx * Ny size
    N = Nx * Ny; % map size

    % Set Tx parameters
    Pt_dBm = 33; % Transmit power = 2000 mW / 33dBm
    freq = 0.8e9; % Frequency = 0.8 GHz / 1.4 GHz / 2.8 GHz

    % ---------------- 导入数据 ----------------
    [DataTable, Tx_tbl] = readData(file_tx, file_rss);
    Power = DataTable.Power;

    % ---------------- 建筑物位置 ----------------
    build = find(sum(Power, 2) <= min(Power(:)) .* width(Power));

    % ------------- 构造发射机结构体 -------------
    Tx_tbl.Properties.VariableNames = {'X','Y','Z'};
    Tx_num = height(Tx_tbl);
    Tx = struct('id',cell(1,Tx_num), 'x',[], 'y',[], 'pos',[], 'Pt_mW', [], 'Pt_dBm',[], 'freq', [], 'gain', [], 'gain_ln',[]);
    for i = 1:Tx_num
        [xi, yi] = raster(Tx_tbl.X(i), Tx_tbl.Y(i), Nx, Ny, unit);
        pos_i = cti(xi, yi, Nx);
        Tx(i).id      = i;
        Tx(i).x       = xi;
        Tx(i).y       = yi;
        Tx(i).pos     = pos_i;

        Tx(i).Pt_dBm  = Pt_dBm;
        Tx(i).Pt_mW   = db_to_ln(Pt_dBm);
        Tx(i).freq    = freq;

        Pr_dBm = Power(:,i);
        gain_dB = Pr_dBm - Pt_dBm;
        Tx(i).gain    = gain_dB;
        Tx(i).gain_ln = db_to_ln(gain_dB);
        
    end

    % -------------- 输出 map 结构 ----------------
    map.name       = mapname;
    map.Nx         = Nx;
    map.Ny         = Ny;
    map.size       = N;
    map.Tx         = Tx;
    map.Tx_num     = Tx_num;
    map.build      = build;

    % Generate phi
    map.phi = phi_ent(map);

end

% Convert dB to linear
function ln = db_to_ln(db)
    ln = 10.^(db./10);
end

% Convert linear to dB
function db = ln_to_db(ln)
    db = 10.*log10(ln);
end

% Convert coordinate to index - 2D
function idx = cti(x, y, Nx)
   idx = (y-1) * Nx + x; 
end

% Convert index to coordinate - 2D
function [x, y] = itc(idx, Nx)
    x = mod(idx - 1, Nx) + 1;
    y = floor((idx - 1) / Nx) + 1;
end

% Rasterization - 2D
function [xi, yi] = raster(x, y, Nx, Ny, unit)
    dx = unit / Nx;  
    dy = unit / Ny;

    xi = floor(x ./ dx) + 1;
    yi = floor(y ./ dy) + 1;

    xi = min(max(xi, 1), Nx);
    yi = min(max(yi, 1), Ny);
end

% Initialize phi_entire without buildings
function [phi] = phi_ent(map)
    N    = map.size;              % map's size
    M    = map.Tx_num;            % Number of Tx
    gains  = [map.Tx.gain_ln];    % N×M
    poses  = [map.Tx.pos];        % 1×M
    
    % phi
    rowIdx = repmat((1:N)', M, 1);
    colIdx = reshape(repmat(poses, N, 1), [], 1 );
    valIdx = gains(:);
    phi = sparse(rowIdx, colIdx, valIdx, N, N);
    phi(map.build, :) = 0;    
end

