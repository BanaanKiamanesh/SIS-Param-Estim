%% Prior Functions
% Calculate the Prior Log Probability for a particular set of parameters.
function p = prior(x)
    A_flat = x(1:end-2);
    beta_ = x(end-2);
    gamma_ = x(end);

    if all(x >= 0)
        % log_p_A =  A_log_prior_symmetric(A_flat);
        log_p_A =  A_log_prior_diagonal(A_flat);
        % log_p_A =  A_log_prior_uniform(A_flat);

        p = log_p_A + beta_log_prior(beta_) + gamma_log_prior(gamma_);
    else
        p = -inf;
    end
end

%% Adjacency Matrix Prior Log Probability Functions
function log_p_A = A_log_prior_symmetric(A_flat)
    A = unflatten(A_flat);
    log_p_A = - norm(A - A', 'fro')^2;;
end

function log_p_A = A_log_prior_diagonal(A_flat)
    A = unflatten(A_flat);
    log_p_A = - norm(A - eye(size(A)), 'fro')^2;
end


% Assuming uniform distribution between 0 and 2 for all components.
function log_p_A = A_log_prior_uniform(A_flat)
    log_p_A = 0;  % for uniform between 0 and 1
    % log_p_A = -(length(A_flat) * log(2));
end

%% Beta Prior Log Probability
function log_p_beta = beta_log_prior(beta_)
    % TODO: Just using a value that would get R = 1.28 with the mean of gamma
    mu = 0.3556;  
    sigma = 0.4;

    log_p_beta = -0.5 * log(2*pi*sigma^2) - ((beta_ - mu)^2) / (2*sigma^2);
end

%% Gamma Prior Log Probability
% Assuming that the Serial Period (1/gamma) is following a gaussian with mean
% on 3.6 days
function log_p_gamma = gamma_log_prior(gamma_)
    mu = 1/3.6;
    sigma = 1/1.6;

    log_p_gamma = -0.5 * log(2*pi*sigma^2) - ((gamma_ - mu)^2) / (2*sigma^2);
end