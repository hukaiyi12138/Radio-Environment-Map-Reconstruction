function psi = psi_mmi(M, N, sample_rate, phi)
    % M: Number of selected SLs (sampling locations)
    % N: Total number of cubes in the ROI
    % phi: Sparse dictionary matrix (size N x N)
    
    % Initialize the set of selected sampling locations (SLs)
    selected_sl = false(1, N);  % Use logical array to track selected locations
    
    % Step 1: Initialize the sensing matrix and its determinant
    current_matrix = sparse(M, N);  % Sparse matrix to store the sensing matrix

    % Step 2: Loop for each sampling location (SL) selection
    samples = round(sample_rate * N);  % Calculate how many samples to select
    for t = 1:samples
        fprintf('Generate psi_mmi: sample_rate = %.2f & process: %d/%d\n', sample_rate, t, samples);

        max_det_value = -inf;
        best_sl = -1;
        
        % Step 3: Try each possible location and calculate its determinant contribution
        for i = 1:N
            if ~selected_sl(i)  % Only consider unselected SLs
                % Construct the candidate sensing matrix by adding the new SL (phi(i,:))
                candidate_matrix = [current_matrix; phi(i, :)];
                
                % Instead of det, we can use a method to incrementally calculate the determinant
                candidate_det = det(candidate_matrix' * candidate_matrix);  % Or use an incremental method
                % For example, using the Cholesky decomposition would be faster
                
                % Check if this is the best SL so far
                if candidate_det > max_det_value
                    max_det_value = candidate_det;
                    best_sl = i;
                end
            end
        end
        
        % Add the best SL to the selected set
        selected_sl(best_sl) = true;  % Mark as selected
        current_matrix = [current_matrix; phi(best_sl, :)];  % Update the sensing matrix with sparse addition
        
        % Optionally, display progress
        fprintf('Selected SL %d, determinant: %.4f\n', best_sl, max_det_value);
        
    end
    
    % Step 4: Construct the final psi matrix
    % Find the indices of selected SLs and construct the final sparse matrix
    col_idx = find(selected_sl);  % Get the indices of selected sampling locations
    row_idx = 1:length(col_idx);  % Row indices corresponding to the selected sampling locations
    
    % Construct the sparse matrix in one go
    psi = sparse(row_idx, col_idx, 1, M, N);  % Final sparse matrix

end
