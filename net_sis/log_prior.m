%% Prior Functions
% Calculate the Prior Log Probability for a particular set of parameters.
function log_p = log_prior(x)
    A_flat = x(1:end-1);
    gamma_ = x(end);

    if all(x >= 0)
        log_p = sum([
            A_log_prior_diagonal(A_flat)
            % A_log_prior_symmetric(A_flat)
            % A_log_prior_uniform(A_flat)
            gamma_log_prior(gamma_)
        ]);
    else
        log_p = -inf;
    end
end

%% Adjacency Matrix Prior Log Probability Functions
% The dynamics of the Network SIS model (in terms of the basic reproduction
% number) are determined by the dominant eigenvalue of A (which encodes the
% transmission rate, typically beta in scalar models). This means that
% there are multiple topologies of the A that result in the same dynamics.
%
% Since the MCMC method is stochastic, it will converge to one of the solutions
% that describes the dynamics but there's is no guarantee that this solution
% will match the correct network topology!
%
% These different functions are an attempt (with limited success) to impart
% some prior knowledge about the actual topology to the algorithm, in the hope
% that it affects the convergence.
function log_p_A = A_log_prior_diagonal(A_flat)
    A = unflatten(A_flat);
    log_p_A = - norm(A - eye(size(A)), 'fro')^2;
end

function log_p_A = A_log_prior_symmetric(A_flat)
    A = unflatten(A_flat);
    log_p_A = - norm(A - A', 'fro')^2;;
end

function log_p_A = A_log_prior_uniform(A_flat)
    log_p_A = 0;  % uniform between 0 and 1
end

%% Gamma Prior Log Probability
% Assuming that the Serial Period (1/gamma) is following a gaussian
function log_p_gamma = gamma_log_prior(gamma_)
    mu = 1/4.1;
    sigma = 0.5;

    log_p_gamma = -0.5 * log(2*pi*sigma^2) - ((gamma_ - mu)^2) / (2*sigma^2);
end