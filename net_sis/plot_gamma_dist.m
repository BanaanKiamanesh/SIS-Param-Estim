function plot_gamma_dist(X, gamma_true, gamma_hat)
    figure;
    title('Simulated \gamma Posterior');
    hold on;
    histogram(X(end,:), 100);
    plot([gamma_true gamma_true], ylim, 'g--', 'LineWidth', 2)
    plot([gamma_hat gamma_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\gamma posterior', '\gamma_{true}', '{\gamma hat}')
end
