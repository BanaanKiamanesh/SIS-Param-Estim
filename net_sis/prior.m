%% Prior Function
% Calculate the Log Probability for a particular set of components, 
% assuming that they are all uniformly distributed between 0 and 1.
function p = prior(x)
    A_flat = x(1:end-1);
    gamma_ = x(end);
    if all(A_flat >= 0 & A_flat <= .5)
        p = A_log_prior(A_flat) + gamma_log_prior(gamma_);
    else
        p = -inf;
    end
end

%% Adjacency Matrix Prior Log Probability
% Assuming uniform distribution between 0 and 0.5 for all components.
function log_p_A = A_log_prior(A_flat)
    log_p_A = -(length(A_flat) * log(0.5));
end

%% Gamma Prior Log Probability
% Assuming that the Serial Period (1/gamma) is following a gaussian with mean
% on 3.6 days
function log_p_gamma = gamma_log_prior(gamma_)
    mu = 1/3.6;
    sigma = 1/1.6;

    log_p_gamma = -0.5 * log(2*pi*sigma^2) - ((gamma_ - mu)^2) / (2*sigma^2);
end
