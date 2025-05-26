clear
close all
clc

%% Parameters

y0 = [0.09 0.15]'; 
tspan = [0 10];
N_nodes = 2;

% True Model Parameters
A_true = [
    1.0 0.3
    0.4 0.8
];
gamma_true = 0.3;

%% Ideal Model
[t, y_true] = sim_net_sis(y0, A_true, gamma_true, tspan);

plot_evolution(t, y_true, 'SIS Model - Population Dynamics')

%% Simulated Noisy Measurements
%Additive Gaussian Noise
sigma = .05;
y_obs = y_true + sigma * randn(size(y_true));

plot_evolution(t, y_obs, 'SIS Model - Noisy Measurements')

%% MCMC Estimation

N_samples = 100000;

% Variance of the Random Walk used for candidate proposal
Var_A = .2 * ones(N_nodes ^ 2, 1);  % For each matrix element
Var_gamma = .2;
Var_c = diag([Var_A; Var_gamma]);  

x = [flatten(eye(N_nodes)); .5];
X = zeros(length(x), N_samples); 
X(:, 1) = x;

% Initial Density Calculation - prior * likelihood
p_current = prior(x) * likelihood(x, t, y_obs, y0, sigma);

f = waitbar(0, 'Running MCMC Estimation...');
for i = 2:N_samples
    % TODO: This generates a lot of candidates with negative elements. 
    % Would it be possible to add an average to `randn` or try a strictly
    % positive distribution? (Would that change any assumption for the method?)
    c = x + Var_c' * randn(size(x));
    if prior(c) == 0  % Avoid the slow likelihood calculation
        X(:, i) = x;  
        continue;
    end

    % Evaluate Density with new candidate parameters
    p_candidate = prior(c) * likelihood(c, t, y_obs, y0, sigma); 

    alpha = min(1, p_candidate / p_current); % Acceptance Probability
    if rand < alpha  % Candidate is accepted (spin the roulette!)
        x = c;
        p_current = p_candidate;
    end
    
    X(:,i) = x;  
    waitbar(i/N_samples, f);
end
close(f)

% %% Minimum Variance Estimate Result
A_hat  = unflatten(mean(X(1:4, :), 2));
gamma_hat = mean(X(end, :));

disp("Results")
disp("A_hat=")
disp(A_hat)
fprintf('Gamma: %0.2f \n', gamma_hat)

plot_X_trajectory(X)
plot_A_dist(X, A_true, A_hat)
plot_gamma_dist(X, gamma_true, gamma_hat)
