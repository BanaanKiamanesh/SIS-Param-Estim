%% Calculate the Log Likelihood
% Calculate the log likelihood for a set of parameters and the observed
% data.
% 
% Assuming a gaussian distribution
function log_L = log_likelihood(x, tspan, I_obs, y0)
    [beta_, kappa, gamma_] = deal(x(1), x(2), x(3));
    try
        [~, y_sim] = simulate_seir(y0, beta_, kappa, gamma_, tspan);
    catch
        disp('Parameters =')
        disp(x)
        error("Could not solve the Likelihood")
    end

    I_sim = y_sim(:, 3);
    residuals = I_obs - I_sim;
    variance = var(residuals);
    n = length(I_obs);

    log_L = -0.5 * n * log(2 * pi * variance) - sum((residuals).^2) / (2 * variance);
end
