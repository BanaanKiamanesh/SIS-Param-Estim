%% Plot the Evolution from the model ODEs
function plot_evolution(t, y, description, leg)
    figure;
    plot(t, y)
    legend(leg)
    xlabel('Time (days)')
    ylabel('Population')
    title(description)
end