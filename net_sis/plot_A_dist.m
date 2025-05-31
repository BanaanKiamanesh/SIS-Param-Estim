%% Plot the Histogram of the posterior 
function plot_A_dist(X, A_true, A_hat)
    figure;
    
    titles = {'A_{11}', 'A_{12}', 'A_{21}', 'A_{22}'};
    
    for j = 1:2
        for i = 1:2
            idx = (j * 2 + i) - 2;
            component = titles{idx};
            subplot(2, 2, idx);
            hold on;
            hist(X(idx, :), 100);
            plot([A_true(i, j), A_true(i, j)], ylim, 'g--', 'LineWidth', 2);
            plot([A_hat(i, j), A_hat(i, j)], ylim, 'r--', 'LineWidth', 2);
            title(['Simulated Posterior of ' component]);
            legend('distribution', component, ['{' component '} hat']);
            hold off;
        end
    end
end