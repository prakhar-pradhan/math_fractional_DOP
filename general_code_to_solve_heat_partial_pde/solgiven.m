% Generalized Solver for Time-Fractional PDEs with L1 Scheme for Example 3.3
alpha_values = [0.2, 0.4, 0.6, 0.8]; % Different alpha values
N_values = 64 * 2.^(0:4); % Mesh sizes: [64, 128, 256, 512, 1024]
T = 1; % Final time
a = pi^2; % Coefficient of u
global u_exact

% Preallocate for results
results = zeros(length(alpha_values), length(N_values), 2); % Errors and rates

for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    errors = zeros(size(N_values));
    
    for idx = 1:length(N_values)
        N = N_values(idx);
        M = N; % Spatial mesh (equal to temporal for simplicity)
        h = pi / M; % Spatial step (domain: x ∈ [0, π])
        tau = T / N; % Temporal step
        
        % Create spatial grid
        x_grid = linspace(0, pi, M+1);
        
        % Define exact solution and source term for y₂
        % u2_exact = @(x, t) t.^2 .* sin(pi*x);

        %% Enter your function here
        f_j_int = @(x, t) (2/gamma(3 - alpha)) * t.^(2 - alpha) .* sin(pi*x) + pi^2 * t.^2 .* sin(pi*x);
        %% 
        
        % L1 Scheme Coefficients
        b = zeros(N, 1);
        for l = 1:N
            b(l) = l^(1 - alpha) - (l-1)^(1 - alpha);
        end
        
        % Initialize Solution - 2D array for time and space
        U2 = zeros(N+1, M+1); % y₂(x,t)
        
        % Initial condition (t=0)
        for j = 1:M+1
            x_j = x_grid(j);
            % This is where you can modify the initial condition for y₂
            %U2(1, j) = 0; % Default: y₂(x,0) = 0 as per the problem
            
            % To set a different initial condition, modify the line above
            %Example: U2(1, j) = sin(x_j); 
            % Example: U2(1, j) = x_j * (pi - x_j);
        end
        
        % Time-Stepping Loop
        for n = 1:N
            t_n = n * tau;
            
            % Space Loop
            for j = 1:M+1
                x_j = x_grid(j);
                
                % Source term at current time and space
                f2_n = f_j_int(x_j, t_n);
                
                % For U2 (y₂)
                if n == 1
                    sum_term_U2 = 0;
                else
                    sum_term_U2 = 0;
                    for l = 1:n-1
                        sum_term_U2 = sum_term_U2 + (b(n-l+1) - b(n-l)) * U2(l+1, j);
                    end
                end
                
                % Update Equation with general coefficient a
                coeff_L1 = b(1) / (gamma(2 - alpha) * tau^alpha);
                denominator = coeff_L1 + a; % Generalized denominator
                
                U2(n+1, j) = (f2_n - sum_term_U2 / (gamma(2 - alpha) * tau^alpha)) / denominator;
            end
        end
        
        % Compute Error at final time for all spatial points
        exact_values_y2 = u_exact(x_grid, T);
        
        % Compute maximum error across all spatial points
        error_y2 = max(abs(U2(N+1, :) - exact_values_y2));
        
        % Use error from y2 solution
        errors(idx) = error_y2;
    end
    
    % Calculate Convergence Rates
    rates = zeros(length(N_values)-1, 1);
    for idx = 1:length(N_values)-1
        rates(idx) = log2(errors(idx)/errors(idx+1));
    end
    
    % Store results
    results(a_idx, :, 1) = errors;
    results(a_idx, 1:end-1, 2) = rates;
end

% Display Results in Table Format
fprintf('\nResults (Maximum Nodal Errors and Rates)\n');
fprintf('Alpha   N       Error         Rate\n');
for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);
    for idx = 1:length(N_values)
        error_val = results(a_idx, idx, 1);
        if idx == 1
            rate_str = '  -  ';
        else
            rate_val = results(a_idx, idx-1, 2);
            rate_str = sprintf('%7.3f', rate_val);
        end
        fprintf('%.1f   %6d  %.3e  %s\n', alpha, N_values(idx), error_val, rate_str);
    end
    fprintf('\n');
end