%% Prior Functions
% Calculate the Prior Log Probability for a particular set of parameters.
function log_p = log_prior(x)
    if all(x >= 0)
        [beta_, kappa, gamma_] = deal(x(1), x(2), x(3));
        log_p = sum([
            beta_log_prior(beta_)
            kappa_log_prior(kappa)
            gamma_log_prior(gamma_)
        ]);
    else
        log_p = -inf;
    end
end

%% Prior Log Probabilities
% Assuming that the prior distributions are gaussian.
% Mean values taken from
%   Chowell G, Miller MA, Viboud C. - Seasonal influenza in the United States,
%   France, and Australia: transmission and prospects for control. 
%   Epidemiol Infect. 2008 Jun;136(6):852-64. doi: 10.1017/S0950268807009144. 
%   Epub 2007 Jul 18. PMID: 17634159; PMCID: PMC2680121. 
%
% Standard Deviations are arbitrary

% Beta Prior Log Probability
function log_p_beta = beta_log_prior(beta_)
    mu = 0.32;
    sigma = 0.5;

    log_p_beta = -0.5 * log(2*pi*sigma^2) - ((beta_ - mu)^2) / (2*sigma^2);
end


% Kappa Prior Log Probability
function log_p_kappa = kappa_log_prior(kappa)
    mu = 1/1.9;
    sigma = 0.5;

    log_p_kappa = -0.5 * log(2*pi*sigma^2) - ((kappa - mu)^2) / (2*sigma^2);
end


% Gamma Prior Log Probability
function log_p_gamma = gamma_log_prior(gamma_)
    mu = 1/4.1;
    sigma = 0.5;

    log_p_gamma = -0.5 * log(2*pi*sigma^2) - ((gamma_ - mu)^2) / (2*sigma^2);
end
