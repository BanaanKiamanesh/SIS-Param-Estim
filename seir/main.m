clear
close all
clc
% Promote ode warning to error, so we can catch it in `likelihood`
warning('error', 'MATLAB:ode45:IntegrationTolNotMet');  

% rng(1296256323)  % For reproducibility

%% Parameters

tspan = [0 200];
y0 = [
    9900 % Susceptible
    90    % Exposed
    10    % Infectious
    0     % Recovered
]';

% True Model Parameters - Trying to be close to a reproduction number R=1.3
beta_true = 0.26;
kappa_true = 0.50;  % Latency Period of 2 days
gamma_true = 0.20;  % Recovery Period of 5 days

fprintf('True Basic Reproduction Number (R_0): %0.4f', beta_true / gamma_true)

%% Ideal Model
[t, y_true] = simulate_seir(y0, beta_true, kappa_true, gamma_true, tspan);

plot_evolution(t, y_true, 'SEIR Model - Population Dynamics', {'S', 'E', 'I', 'R'})

%% Simulated Noisy Measurements
% Additive Gaussian Noise
sigma = 25;
I_true = y_true(:, 3);
I_obs = I_true + sigma * randn(size(I_true));

plot_evolution(t, [I_true I_obs], 'Noisy Measurements', {'True Infections', 'Observed Infections'})

%% MCMC Estimation

% N_samples = 10000;
N_samples = 100000;

% Initial Estimate
x = [0.31 0.53 0.24]';  % Prior Mean

% Variance of the Random Walk used for candidate proposal
Var_c = 0.0016 * eye(length(x));

% Memory Allocation
X = zeros(length(x), N_samples); 
X(:, 1) = x;

Alpha = zeros(1, N_samples);

% Initial Density Calculation - prior * likelihood
p_current = log_prior(x) + log_likelihood(x, t, I_obs, y0);

f = waitbar(0);
accepted = 0;
for i = 2:N_samples
    c = x + Var_c * randn(size(x));
    if log_prior(c) == -inf  % Avoid slow likelihood calculation on 0 probability
        X(:, i) = x;  
        continue;
    end

    % Evaluate Density with new candidate parameters
    p_candidate = log_prior(c) + log_likelihood(c, t, I_obs, y0);

    alpha = min(1, exp(p_candidate - p_current)); % Acceptance Probability
    Alpha(i) = alpha;
    if rand < alpha  % Candidate is accepted (spin the roulette!)
        x = c;
        p_current = p_candidate;
        accepted = accepted + 1;
    end
    
    X(:,i) = x; 
    waitbar(i/N_samples, f, sprintf('MCMC - Acp. Rate: %0.2f', accepted / i));
end
close(f)

% Has to be beetwen 30-40% on average
acceptance_rate_percent = (accepted / N_samples) * 100; 

%% Minimum Variance Estimate Results
burn_in = floor(0.15 * size(X, 2));  % burn-in of 15%
X_post = X(:, burn_in + 1:end);       

beta_hat = mean(X_post(1, :));
kappa_hat = mean(X_post(2, :));
gamma_hat = mean(X_post(3, :));

R_hat = beta_hat / gamma_hat;

std_beta_hat = std(X_post(1, :), 0, 2);       
std_kappa_hat = std(X_post(2, :), 0, 2);       
std_gamma_hat = std(X_post(3, :), 0, 2);       

disp("Results")
disp(['Acceptance rate: ', num2str(acceptance_rate_percent), '%']);
disp("---------------")
fprintf('beta_hat: %0.4f \n', beta_hat)
fprintf('kappa_hat: %0.4f \n', kappa_hat)
fprintf('gamma_hat: %0.4f \n', gamma_hat)
fprintf('Estimated Effective Reproduction Number (R): %0.4f \n', R_hat)
disp("---------------")
fprintf('std_beta_hat: %0.4f \n', std_beta_hat)
fprintf('std_kappa_hat: %0.4f \n', std_kappa_hat)
fprintf('std_gamma_hat: %0.4f \n', std_gamma_hat)

plot_X_trajectory(X, burn_in)
plot_distributions(...
    X_post, ...
    beta_true, ...
    beta_hat, ...
    kappa_true, ...
    kappa_hat, ...
    gamma_true, ...
    gamma_hat ...
)
