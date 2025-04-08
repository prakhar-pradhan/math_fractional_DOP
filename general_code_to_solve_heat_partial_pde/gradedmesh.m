function [U, t, x] = solve_fractional_diffusion(alpha, T, X_val, M, N, r)
    % Generate time mesh
    t = zeros(N+1,1);
    for k=0:N
        t(k+1) = T * (k/N)^r;
    end
    % Generate spatial mesh
    x = linspace(0, X_val, M+1);
    h = x(2) - x(1);
    % Initialize solution
    U = zeros(N+1, M+1);
    U(1,:) = sin(x);
    % Precompute Gamma(2 - alpha)
    Gamma_2_minus_alpha = gamma(2 - alpha);
    
    % Compute spatial operator matrix A for interior points (m=2 to M)
    A = zeros(M-1, M-1);
    for i=1:M-1
        A(i,i) = 2/h^2 + (1 + x(i+1));
        if i > 1
            A(i,i-1) = -1/h^2;
        end
        if i < M-1
            A(i,i+1) = -1/h^2;
        end
    end
    
    % Time stepping
    for n=1:N
        j = n+1;
        tau_n = t(j) - t(j-1);
        % Compute b_k for k=0 to n-1
        b = zeros(n,1);
        for k=0:n-1
            tau_kp1 = t(k+2) - t(k+1);
            if tau_kp1 == 0
                continue;
            end
            term1 = (t(j) - t(k+1))^(1-alpha);
            term2 = (t(j) - t(k+2))^(1-alpha);
            b(k+1) = 1 / Gamma_2_minus_alpha / tau_kp1 * (term1 - term2);
        end
        % Compute a_vec for coefficients of past terms
        a_vec = zeros(n+1,1); % For U(1) to U(n+1)
        a_vec(1) = -b(1); % For l=0, coefficient of U(1,m)
        for l=2:n+1
            if l <= n
                a_vec(l) = b(l-1) - b(l); % For l=1 to n, coefficient of U(l+1,m)
            else
                % For l=n+1, coefficient of U(n+1,m) is b(n), from earlier
                a_vec(n+1) = b(n); % Since for l=n, coefficient is b(n-1) - b(n), but in our sum, it's b(n-1)
            end
        end
        % Compute past part P_j_int
        U_hist = U(1:n,:);
        P_j_full = a_vec(1:n)' * U_hist(:,2:M); % Sum up to U(n,:), exclude U(n+1,:)
        % Compute f_j_int
        %f_j_int = (x(2:M).(pi - x(2:M)).(1 + t(j)^4) + t(j)^2)';
        %% Enter your function here
        f_j_int = (x(2:M) .* (pi - x(2:M)) .* (1 + t(j)^4) + t(j)^2)';
        %%

        
        % Coefficient of U(n+1,m) in D_N^alpha is b(n), which is a_vec(n+1)
        c_j = a_vec(n+1); % Coefficient of U(n+1,m)
        % Solve linear system
        matrix = c_j * eye(M-1) + A;
        rhs = f_j_int - P_j_full';
        U_j_int = matrix \ rhs;
        % Update solution
        U(j,2:M) = U_j_int';
        U(j,1) = 0;
        U(j,M+1) = 0;
    end
end

function D = compute_two_mesh_difference(U_coarse, U_fine, t_coarse, t_fine, x_coarse, x_fine)
    [Nc, Mc] = size(U_coarse);
    [Nf, Mf] = size(U_fine);
    if Mf ~= 2*Mc - 1 || Nf ~= 2*Nc - 1
        error('Fine mesh must be twice the coarse mesh in both dimensions');
    end
    max_diff = 0;
    for i=1:Nc
        for j=1:Mc
            % Interpolate fine solution to coarse grid points
            t_c = t_coarse(i);
            x_c = x_coarse(j);
            % Find nearest points in fine grid
            t_f_idx = find(abs(t_fine - t_c) == min(abs(t_fine - t_c)), 1);
            x_f_idx = find(abs(x_fine - x_c) == min(abs(x_fine - x_c)), 1);
            diff_val = abs(U_coarse(i,j) - U_fine(t_f_idx, x_f_idx));
            if diff_val > max_diff
                max_diff = diff_val;
            end
        end
    end
    D = max_diff;
end

% Main script for both graph and table
clear all;
alpha = 0.6;
T = 1;
X_val = pi;
M_plot = 100;
N_plot = 100;
r_plot = (2 - alpha)/alpha;

% Plot the graph
[U_plot, t_plot, x_plot] = solve_fractional_diffusion(alpha, T, X_val, M_plot, N_plot, r_plot);
figure(1);
[ X_grid, Y_grid ] = meshgrid(x_plot, t_plot');
surf(X_grid, Y_grid, U_plot');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Solution u(x,t) for \alpha = 0.6');

% Reproduce Table
r_values = [1, (2 - alpha)/(2*alpha), (2 - alpha)/alpha, 2*(2 - alpha)/alpha];
 M_list = [64; 128; 256; 512; 1024; 2048];
%M_list = [4];
fprintf('Table 7 for alpha = 0.6:\n');
fprintf('r\tN=M=64\tN=M=128\tN=M=256\tN=M=512\tN=M=1024\n');
for k=1:length(r_values)
    r = r_values(k);
    fprintf('r=%g\n', r);
    D_list = zeros(1, length(M_list)-1);
    for j=1:length(M_list)-1
        M_coarse = M_list(j);
        N_coarse = M_coarse;
        M_fine = M_list(j+1);
        N_fine = M_fine;
        [U_coarse, t_coarse, x_coarse] = solve_fractional_diffusion(alpha, T, X_val, M_coarse, N_coarse, r);
        [U_fine, t_fine, x_fine] = solve_fractional_diffusion(alpha, T, X_val, M_fine, N_fine, r);
        D_list(j) = compute_two_mesh_difference(U_coarse, U_fine, t_coarse, t_fine, x_coarse, x_fine);
    end
    % Compute orders for N=64 to 512
    p_list = zeros(1, length(M_list)-2);
    for i=1:length(M_list)-2
        p_list(i) = log2(D_list(i) / D_list(i+1));
    end
    % Print results
    fprintf('D: ');
    for i=1:length(D_list)
        fprintf('%e\t', D_list(i));
    end
    fprintf('\n');
    fprintf('p: ');
    for i=1:length(p_list)
        fprintf('%f\t', p_list(i));
    end
    fprintf('\n\n');
end