function plot_beta_gamma_dist(X, beta_true, beta_hat, gamma_true, gamma_hat)
    figure;

    subplot(1, 2, 1);
    title('Simulated \beta Posterior');
    hold on;
    hist(X(end-1,:), 100);
    plot([beta_true beta_true], ylim, 'g--', 'LineWidth', 2)
    plot([beta_hat beta_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\beta posterior', '\beta_{true}', '{\beta hat}')

    subplot(1, 2, 2);
    title('Simulated \gamma Posterior');
    hold on;
    hist(X(end,:), 100);
    plot([gamma_true gamma_true], ylim, 'g--', 'LineWidth', 2)
    plot([gamma_hat gamma_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\gamma posterior', '\gamma_{true}', '{\gamma hat}')
end
