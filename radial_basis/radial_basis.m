% MATLAB code to solve a fractional differential equation using RBFs
% Equation: D^0.8 y(x) = (Γ(3)/Γ(2.2)) * x^1.2, with y(0) = 0
% Exact solution: y(x) = x^2

%% RBF Kernel Function Definition
function phi = rbf(type, r, beta, nu)
    % Returns value of the RBF based on type and input radius r
    switch type
        case 'gaussian'
            phi = exp(-r^2 / 2);
        case 'multiquadric'
            phi = (1 + r^2)^(beta/2);
        case 'matern'
            c = 10; % Scaling factor to stabilize the Bessel function
            phi = (c*r)^nu * besselk(nu, c*r); 
        case 'powers'
            phi = r^beta;
        case 'thin_plate'
            phi = r^(2*beta) * log(r + 1e-20); % Add small constant to avoid log(0)
        otherwise
            error('Unknown RBF type');
    end
end

%% Modified Basis Function (centered at x_tilde(j))
function p = phi_j(x, j, type, beta, nu, x_tilde, c_factor)
    r = abs((x - x_tilde(j))/(c_factor)); % Scaled distance
    if strcmp(type, 'matern') && r < 1e-20
        r = 1e-20; % Avoid singularity in Bessel function at r = 0
    end
    p = rbf(type, r, beta, nu); % Evaluate RBF
end

%% Grünwald–Letnikov Approximation for Fractional Derivative
function df = gl_fd(f, alpha, x, h)
    m = floor(x / h); % Number of steps in the past
    df = 0;
    for k = 0:m
        % Compute binomial coefficient for fractional order
        binom = 1;
        for j = 1:k
            binom = binom * (alpha - j + 1) / j;
        end
        t = x - k * h;
        if t < 0
            continue % Outside domain
        end
        df = df + (-1)^k * binom * f(t); % Weighted sum
    end
    df = df / h^alpha; % Scale by h^alpha
end

%% MAIN SCRIPT
n = 50; % Number of collocation points
c_factor = 1; % Scaling factor for radial distance
x = linspace(0,1,n); % Collocation points in [0, 1]
x_tilde = x; % Centers for RBFs

alpha = 0.8; % Order of fractional derivative
h = 1e-4; % Step size for GL method
reg_param = 1e-12; % Tikhonov regularization for numerical stability

% List of RBF types and corresponding parameters
types = {'gaussian', 'multiquadric', 'matern', 'powers', 'thin_plate'};
betas = [1, 1, 1.5, 3, 2]; % β values
nus = [0, 0, 1.5, 0, 0];   % ν values (used for matern only)

% Preallocate output cell array
y_approx_all = cell(1,5);

% Loop through each RBF type
for type_idx = 1:5
    type = types{type_idx};
    beta_val = betas(type_idx);
    nu_val = nus(type_idx);
    
    % Build matrix C of size n × n
    C = zeros(n,n);
    
    % First row enforces initial condition y(0) = 0
    for j = 1:n
        C(1,j) = phi_j(x(1), j, type, beta_val, nu_val, x_tilde, c_factor);
    end
    
    % Populate remaining rows with fractional derivative evaluations
    for i = 2:n
        for j = 1:n
            phi_j_func = @(t) phi_j(t, j, type, beta_val, nu_val, x_tilde, c_factor);
            C(i,j) = gl_fd(phi_j_func, alpha, x(i), h); % Apply GL derivative
        end
    end
    
    % Add regularization term for numerical stability
    C = C + reg_param * eye(n);
    
    % Right-hand side of equation
    G = zeros(n, 1);
    G(2:n) = (gamma(3)/gamma(2.2)) * x(2:n).^1.2; % f(x) for x > 0

    % Solve the linear system C * lambda = G
    try
        lambda = C \ G;
    catch ME
        warning(['Failed to solve for ' type ' RBF: ' ME.message]);
        y_approx_all{type_idx} = NaN(size(linspace(0,1,200)));
        continue;
    end
    
    % Evaluate solution y(x) = ∑ λ_j * φ_j(x)
    x_plot = linspace(0,1,200);
    y_approx = zeros(size(x_plot));
    for k = 1:length(x_plot)
        for j = 1:n
            y_approx(k) = y_approx(k) + lambda(j) * phi_j(x_plot(k), j, type, beta_val, nu_val, x_tilde, c_factor);
        end
    end
    y_approx_all{type_idx} = y_approx;
    
    % Compare to exact solution y(x) = x^2
    y_exact = x_plot.^2;
    if all(~isnan(y_approx))
        error_inf = max(abs(y_exact - y_approx));
        fprintf('L-infinity error for %-15s: %.6e\n', type, error_inf);
    end
end

%% Plotting the results
figure;
x_plot = linspace(0,1,200);
y_exact = x_plot.^2; % Exact solution
plot(x_plot, y_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Exact');
hold on;

colors = ['r', 'g', 'c', 'm', 'k']; % Colors for each RBF
line_styles = ['--', '-.', ':', '-', '--']; % Line styles

for type_idx = 1:5
    if ~isnan(y_approx_all{type_idx}(1))
        type = types{type_idx};
        plot(x_plot, y_approx_all{type_idx}, [colors(type_idx) line_styles(type_idx)],...
            'LineWidth', 2, 'DisplayName', [type ' Approximation']);
    end
end

hold off;
legend('Location', 'northwest');
title('Solution using different RBFs for D^{0.8} y(x) = (Γ(3)/Γ(2.2)) x^{1.2}');
xlabel('x');
ylabel('y(x)');
grid on;
