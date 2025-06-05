%% Ideal evolution of the Network SIS Model for a set of parameters.
% y0     - Initial Conditions [n x 1]
% A      - Contact Network Adjacency Matrix [n x n]
% betta_ - Infection Rate
% gamma_ - Recovery Rate
% tspan  - Simulation Time
function [t, y_sim] = sim_net_sis(y0, A, gamma_, tspan)
    [t, y_sim] = ode45(@system_ode, tspan, y0);

    function dxdt = system_ode(t, x)
        dxdt = diag(ones(size(x)) - x) * A * x - gamma_ * x;
    end
end
