%% Plot the Evolution from the model ODEs
function plot_evolution(t, y, description)
    figure;
    plot(t, y)
    legend('y_1', 'y_2')
    xlabel('Time (days)')
    ylabel('Population')
    title(description)
end