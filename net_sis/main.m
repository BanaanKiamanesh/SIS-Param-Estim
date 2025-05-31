clear
close all
clc

%% Parameters

dt = 0.1;
tspan = 0:dt:10;
y0 = [0.09 0.15]';
N_nodes = 2;

% True Model Parameters
A_true = [
    1.5 0.3
    0.4 0.8
];
gamma_true = 0.5;

%% Ideal Model
[t, y_true] = sim_net_sis(y0, A_true, gamma_true, tspan);

plot_evolution(t, y_true, 'SIS Model - Population Dynamics')

%% Simulated Noisy Measurements
% Additive Gaussian Noise
sigma = 0.015;
y_obs = y_true + sigma * randn(size(y_true));

plot_evolution(t, y_obs, 'SIS Model - Noisy Measurements')

%% MCMC Estimation

N_samples = 100000;

% Variance of the Random Walk used for candidate proposal
Var_A = .018 * ones(N_nodes ^ 2, 1);  % For each matrix element
Var_gamma = .018;
Var_c = diag([Var_A; Var_gamma]);  

x = [flatten(eye(N_nodes)); .5];
X = zeros(length(x), N_samples); 
X(:, 1) = x;

% Initial Density Calculation - prior * likelihood
p_current = prior(x) + likelihood(x, t, y_obs, y0, sigma);

f = waitbar(0, 'Running MCMC Estimation...');
accepted = 0;
for i = 2:N_samples
    c = x + Var_c' * randn(size(x));
    if prior(c) == 0  % Avoid the slow likelihood calculation
        X(:, i) = x;  
        continue;
    end

    % Evaluate Density with new candidate parameters
    p_candidate = prior(c) + likelihood(c, t, y_obs, y0, sigma); 

    alpha = min(1, exp(p_candidate - p_current)); % Acceptance Probability
    if rand < alpha  % Candidate is accepted (spin the roulette!)
        x = c;
        p_current = p_candidate;
        accepted = accepted + 1;
    end
    
    X(:,i) = x; 
    waitbar(i/N_samples, f);
end
close(f)

% Has to be beetwen 30-40% on average
acceptance_rate_percent = (accepted / N_samples) * 100; 

%% Minimum Variance Estimate Results
burn_in = floor(0.15 * size(X, 2));  % burn-in of 15%
X_post = X(:, burn_in+1:end);       

A_hat     = unflatten(mean(X_post(1:4, :), 2));
gamma_hat = mean(X_post(end, :));

std_A_hat = unflatten(std(X_post(1:4, :), 0, 2));  
std_gamma  = std(X_post(end, :), 0, 2);       

disp("Results")
disp(['Acceptance rate: ', num2str(acceptance_rate_percent), '%']);
disp("---------------")
disp("A_hat=")
disp(A_hat)
fprintf('Gamma: %0.2f \n', gamma_hat)
disp("---------------")
disp("std_A_hat=")
disp(std_A_hat)
fprintf('std_gamma: %0.2f \n', std_gamma)

plot_X_trajectory(X)
plot_A_dist(X_post, A_true, A_hat)
plot_gamma_dist(X_post, gamma_true, gamma_hat)
