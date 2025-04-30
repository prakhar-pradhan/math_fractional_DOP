% Generalized Solver for Time-Fractional PDEs with L1 Scheme for Example 3.3
% Case when exact solution is given
% Define fractional orders to test
alpha_values = [0.2, 0.4, 0.6, 0.8];

% Define spatial/temporal resolutions (refinement levels)
N_values = 64 * 2.^(0:4); % N = [64, 128, 256, 512, 1024]

% Final simulation time
T = 1;

% Coefficient a from the PDE (e.g., from -a*u_xx), often pi^2
a = pi^2;

% Assume u_exact(x,t) is defined elsewhere globally
global u_exact

% Initialize a 3D array to store errors and convergence rates
% Dimensions: (alpha_idx, N_idx, [error, rate])
results = zeros(length(alpha_values), length(N_values), 2);

% Loop over different alpha (fractional derivative order) values
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    errors = zeros(size(N_values)); % To store errors for each mesh size
    
    % Loop over different mesh sizes
    for idx = 1:length(N_values)
        N = N_values(idx);  % Number of time steps
        M = N;              % Number of spatial steps (same as time for simplicity)
        h = pi / M;         % Spatial step size over [0, π]
        tau = T / N;        % Temporal step size
        
        % Create spatial grid over domain [0, π]
        x_grid = linspace(0, pi, M+1);
        
        % Define the source term f(x,t) for y2 equation
        % Example: for exact solution u2 = t^2 * sin(pi*x), the source is:
        f_j_int = @(x, t) (2/gamma(3 - alpha)) * t.^(2 - alpha) .* sin(pi*x) ...
                        + pi^2 * t.^2 .* sin(pi*x);
        
        % Precompute L1 weights (convolution coefficients)
        b = zeros(N, 1);
        for l = 1:N
            b(l) = l^(1 - alpha) - (l-1)^(1 - alpha);
        end
        
        % Initialize numerical solution U2 over (time x space)
        U2 = zeros(N+1, M+1);
        
        % Set initial condition at t = 0 (here, assumed zero)
        for j = 1:M+1
            x_j = x_grid(j);
            % Modify below to apply a different initial condition
            U2(1, j) = 0;
        end
        
        % Time stepping loop
        for n = 1:N
            t_n = n * tau;
            
            % Loop over all spatial points
            for j = 1:M+1
                x_j = x_grid(j);
                
                % Evaluate source term at (x_j, t_n)
                f2_n = f_j_int(x_j, t_n);
                
                % Compute convolution sum for Caputo derivative
                if n == 1
                    sum_term_U2 = 0;
                else
                    sum_term_U2 = 0;
                    for l = 1:n-1
                        sum_term_U2 = sum_term_U2 + ...
                            (b(n-l+1) - b(n-l)) * U2(l+1, j);
                    end
                end
                
                % Coefficient from L1 scheme for time-fractional derivative
                coeff_L1 = b(1) / (gamma(2 - alpha) * tau^alpha);
                
                % Solve for U2 at next time step using the implicit formula
                denominator = coeff_L1 + a;
                U2(n+1, j) = (f2_n - sum_term_U2 / ...
                             (gamma(2 - alpha) * tau^alpha)) / denominator;
            end
        end
        
        % Evaluate exact solution at final time for error calculation
        exact_values_y2 = u_exact(x_grid, T);
        
        % Compute maximum nodal error across space at final time
        error_y2 = max(abs(U2(N+1, :) - exact_values_y2));
        errors(idx) = error_y2;
    end
    
    % Compute convergence rates using log2 ratio of successive errors
    rates = zeros(length(N_values)-1, 1);
    for idx = 1:length(N_values)-1
        rates(idx) = log2(errors(idx)/errors(idx+1));
    end
    
    % Store errors and rates in the results matrix
    results(a_idx, :, 1) = errors;
    results(a_idx, 1:end-1, 2) = rates;
end

% Display tabulated results for each alpha and mesh size
fprintf('\nResults (Maximum Nodal Errors and Rates)\n');
fprintf('Alpha   N       Error         Rate\n');
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    for idx = 1:length(N_values)
        error_val = results(a_idx, idx, 1);
        if idx == 1
            rate_str = '  -  '; % No rate for first N
        else
            rate_val = results(a_idx, idx-1, 2);
            rate_str = sprintf('%7.3f', rate_val);
        end
        fprintf('%.1f   %6d  %.3e  %s\n', alpha, N_values(idx), error_val, rate_str);
    end
    fprintf('\n');
end
