%% Ideal evolution of the Network SIS Model for a set of parameters.
% y0     - Initial Conditions [n x 1]
% N      - Population Size
% beta_  - Infection Rate
% gamma_ - Recovery Rate
% kappa  - Latency Rate
% tspan  - Simulation Time

function [t, y_sim] = simulate_seir(y0, beta_, kappa, gamma_, tspan)
    [t, y_sim] = ode45(@system_ode, tspan, y0);

    function dydt = system_ode(t, y)
        [S, E, I, P] = deal(y(1), y(2), y(3), y(4));
        N = sum(y);
        dydt = [
            -beta_ * S * (I / N)
            beta_ * S * (I / N) - kappa * E
            kappa * E - (gamma_) * I
            gamma_ * I
        ];
    end
end

