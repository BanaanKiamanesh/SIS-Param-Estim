%% Prior Functions
% Calculate the Prior Log Probability for a particular set of parameters.
function log_p = log_prior(x)
    if all(x >= 0 & x <= 1)
        log_p = 0;
        % [beta_, kappa, gamma_] = deal(x(1), x(2), x(3));
        % log_p = sum([
        %     beta_log_prior(beta_)
        %     kappa_log_prior(kappa)
        %     gamma_log_prior(gamma_)
        % ]);
    else
        log_p = -inf;
    end
end

%% Beta Prior Log Probability
function log_p_beta = beta_log_prior(beta_)
    % TODO: Adjust values
    mu = 2;
    sigma = 0.4;

    log_p_beta = -0.5 * log(2*pi*sigma^2) - ((beta_ - mu)^2) / (2*sigma^2);
end

%% Kappa Prior Log Probability
function log_p_kappa = kappa_log_prior(kappa)
    mu = 1/1.9;
    sigma = 0.01;

    log_p_kappa = -0.5 * log(2*pi*sigma^2) - ((kappa - mu)^2) / (2*sigma^2);
end


%% Gamma Prior Log Probability
% Assuming that the Serial Period (1/gamma) is following a gaussian
function log_p_gamma = gamma_log_prior(gamma_)
    mu = 1/4.1;
    sigma = 1/1.6;  % TODO: Adjust

    log_p_gamma = -0.5 * log(2*pi*sigma^2) - ((gamma_ - mu)^2) / (2*sigma^2);
end
