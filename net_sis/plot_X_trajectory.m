%% MCMC Trajectory Plot
function plot_X_trajectory(X)
    figure;
    title('MCMC Trajectory')

    A_11 = X(1, :);
    gamma_ = X(end, :);

    % Plot Trajectory
    hold on;
    plot(A_11, gamma_, 'o--', A_11(1), gamma_(1), 'r*')

    % Fit a 2D Gaussian to the data and draw ellipsis highlighting 1, 2 and
    % 3 std. deviation.
    mu = mean([A_11', gamma_']);
    Sigma = cov([A_11', gamma_']);

    plot_gaussian_ellipse(mu, Sigma, 1, 'g');
    plot_gaussian_ellipse(mu, Sigma, 2, 'y');
    plot_gaussian_ellipse(mu, Sigma, 3, 'r');

    plot(mu(1), mu(2), 'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
end

function plot_gaussian_ellipse(mu, Sigma, nsig, color)
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];
    [V, D] = eig(Sigma);
    ellipse = nsig * V * sqrt(D) * circle;
    ellipse = ellipse + mu';
    plot(ellipse(1,:), ellipse(2,:), color, 'LineWidth', 0.5);
end