clear
close all
clc
% Promote ode warning to error, so we can catch it in `likelihood`
warning('error', 'MATLAB:ode45:IntegrationTolNotMet');  

rng(1296256323)  % For report reproducibility

%% Parameters

dt = 1;
tspan = 0:dt:120;
y0 = [0.01 0.018]';
N_nodes = 2;

% True Model Parameters - Trying to be close to a reproduction number R=1.28,
% inspired in the common cold (influenza)).
A_true = [
  0.7  0.3
  0.4  0.7
];
beta_true = 0.26;
gamma_true = 0.2;  % Recovery Period of 5 days

fprintf('True Effective Reproduction Number (R): %0.4f', calculate_R(A_true, beta_true, gamma_true))

%% Ideal Model
[t, y_true] = sim_net_sis(y0, A_true, beta_true, gamma_true, tspan);

plot_evolution(t, y_true, 'SIS Model - Population Dynamics')

%% Simulated Noisy Measurements
% Additive Gaussian Noise
sigma = 0.015;
y_obs = y_true + sigma * randn(size(y_true));
y_obs = abs(y_obs);  % Ensure only positive values.

plot_evolution(t, y_obs, 'SIS Model - Noisy Measurements')


%% MCMC Estimation

% N_samples = 10000;
N_samples = 100000;
% N_samples = 1000000;

% Variance of the Random Walk used for candidate proposal
Var_A = 0.01 * ones(N_nodes^2, 1);
Var_beta = 0.001;
Var_gamma = 0.001;
Var_c = diag([Var_A; Var_beta; Var_gamma]);  

% Initial Estimate
x = [(1/N_nodes) * ones(N_nodes^2, 1); 0.25; 0.25];

% Memory Allocation
X = zeros(length(x), N_samples); 
X(:, 1) = x;

% Initial Density Calculation - prior * likelihood
p_current = prior(x) + likelihood(x, t, y_obs, y0, sigma);

f = waitbar(0, 'Running MCMC Estimation...');
accepted = 0;
for i = 2:N_samples
    c = x + Var_c' * randn(size(x));
    if prior(c) == -inf  % Avoid slow likelihood calculation on 0 probability
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
X_post = X(:, burn_in + 1:end);       

A_hat = unflatten(mean(X_post(1:N_nodes^2, :), 2));
beta_hat = mean(X_post(end-1, :));
gamma_hat = mean(X_post(end, :));

R_hat = calculate_R(A_hat, beta_hat, gamma_hat);

std_A_hat = unflatten(std(X_post(1:N_nodes^2, :), 0, 2));  
std_beta_hat = std(X_post(end-1, :), 0, 2);       
std_gamma_hat = std(X_post(end, :), 0, 2);       

disp("Results")
disp(['Acceptance rate: ', num2str(acceptance_rate_percent), '%']);
disp("---------------")
disp("A_hat=")
disp(A_hat)
fprintf('beta_hat: %0.4f \n', beta_hat)
fprintf('gamma_hat: %0.4f \n', gamma_hat)
fprintf('Estimated Effective Reproduction Number (R): %0.4f \n', R_hat)
disp("---------------")
disp("std_A_hat=")
disp(std_A_hat)
fprintf('std_beta_hat: %0.4f \n', std_beta_hat)
fprintf('std_gamma_hat: %0.4f \n', std_gamma_hat)

plot_X_trajectory(X, burn_in)
plot_A_dist(X_post, A_true, A_hat)
plot_beta_gamma_dist(X_post, beta_true, beta_hat, gamma_true, gamma_hat)

[t, y_hat] = sim_net_sis(y0, A_hat, beta_hat, gamma_hat, tspan);
plot_evolution(t, y_true - y_hat, 'Evolution Estimation Error')