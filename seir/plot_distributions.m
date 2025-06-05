function plot_distributions(X, beta_true, beta_hat, kappa_true, kappa_hat, gamma_true, gamma_hat)
    figure;
    subplot(2, 2, 1)
    title('Simulated \beta Posterior');
    hold on;
    histogram(X(1,:), 100);
    plot([beta_true beta_true], ylim, 'g--', 'LineWidth', 2)
    plot([beta_hat beta_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\beta posterior', '\beta_{true}', '{\beta hat}')

    subplot(2, 2, 2)
    title('Simulated \kappa Posterior');
    hold on;
    histogram(X(2,:), 100);
    plot([kappa_true kappa_true], ylim, 'g--', 'LineWidth', 2)
    plot([kappa_hat kappa_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\kappa posterior', '\kappa_{true}', '{\kappa hat}')

    subplot(2, 2, 3)
    title('Simulated \gamma Posterior');
    hold on;
    histogram(X(3,:), 100);
    plot([gamma_true gamma_true], ylim, 'g--', 'LineWidth', 2)
    plot([gamma_hat gamma_hat], ylim, 'r--', 'LineWidth', 2)
    hold off;
    legend('\gamma posterior', '\gamma_{true}', '{\gamma hat}')
end
