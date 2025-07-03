%% OMP algorithm
% y: Observation vector
% A: Sensing matrix
% t: Number of iterations
function [theta, rn_norm, converge_point] = omp(y, A, t)
    % Set the convergence tolerance
    tol = 1e-6;

    % Pre-allocate storage space
    rn_norm = zeros(1, t);  % Store the norm of the residual r_n at each iteration
    converge_point = 1; % Initialize the convergence point

    [y_rows, y_columns] = size(y);
    if y_rows < y_columns
        y = y';     % Ensure y is a column vector
    end
    [m, n] = size(A);  % A is m x n

    theta = zeros(n, 1);  % Store the estimated sparse signal (theta)
    theta_ls_current = theta;
    At = zeros(m, t);     % Store selected columns of A during iterations
    Pos_theta = zeros(1, t);  % Store the indices of selected columns of A during iterations
    r_n = y;              % Initialize the residual to y

    for i = 1:t
        product = A' * r_n;  % Compute the inner product between each column of A and the residual
        [~, pos] = max(abs(product));  % Find the column with the highest absolute inner product / the most correlated with the residual
        At(:, i) = A(:, pos);         % Store this column
        Pos_theta(i) = pos;           % Store the index of this column
        A(:, pos) = zeros(m, 1);      % Set this column of A to zero (orthogonalize)

        % Use lsqnonneg to ensure that theta is non-negative
        theta_ls = lsqnonneg(At(:, 1:i), y);% Solve the least squares problem ensuring theta is non-negative
        r_n = y - At(:, 1:i) * theta_ls;  % Update the residual

        % Compute the norm of the residual
        rn_norm(i) = norm(r_n);

        % Check for convergence
        converge_point = i;
        if i >= 2 && abs(rn_norm(i) - rn_norm(i-1)) <= tol
            converge_point = i-1;
            theta_ls = theta_ls_current;
            break;
        end
        theta_ls_current = theta_ls;
    end

    % Recover the final theta
    for j = 1:converge_point
        theta(Pos_theta(j)) = theta_ls(j);  % Fill the positions of theta with the corresponding values from theta_ls
    end

end
