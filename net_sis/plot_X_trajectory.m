%% MCMC Trajectory Plot
function plot_X_trajectory(X)
    figure;
    title('MCMC Trajectory')
    % Restrict to 1st component of A and gamma
    plot(X(1, :), X(end, :), 'o--', X(1,1), X(end,1), 'r*')
end