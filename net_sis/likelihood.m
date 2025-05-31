%% Calculate the Likelihood
% Calculate the log likelihood for a set of parameters and the observed
% data.
% 
% Assuming a gaussian distribution
function L = likelihood(theta, tspan, y_obs, y0, sigma)
    A = unflatten(theta(1:end - 1));
    gamma_ = theta(end);

    try
        [~, y_sim] = sim_net_sis(y0, A, gamma_, tspan);
    catch
        disp('A=')
        disp(A)
        disp('gamma=')
        disp(gamma_)
        error("Could not solve the Likelihood")
    end

    residual = y_obs - y_sim;
    n = numel(y_obs);
    
    L = -0.5 * n * log(2 * pi * sigma^2) - sum(residual(:).^2) / (2 * sigma^2);
end
