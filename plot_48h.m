function plot_48h(tspan, y, idx, var_name, y_patch_max, color)

    t_range = (tspan >= 656) & (tspan <= 704);     %Extract 48-hour window
    t_48h   = tspan(t_range) - 656; 
    var_48h = y(t_range, idx);

    dark_periods = [12, 24; 36, 48];
    labels = {'Dark Period 1 (7pm - 7am)', 'Dark Period 2 (7pm - 7am)'};

    figure;
    ax = axes;
    hold on;

    for i = 1:size(dark_periods, 1)
        start_t = dark_periods(i, 1);
        end_t   = dark_periods(i, 2);
        x_patch = [start_t, end_t, end_t, start_t];
        y_patch = [0, 0, y_patch_max, y_patch_max];
        patch(ax, x_patch, y_patch, [0.5 0.5 0.5],'FaceAlpha', 0.3, 'EdgeColor', 'black', 'DisplayName', labels{i});
    end

    plot(ax, t_48h, var_48h, color, 'LineWidth', 2, 'DisplayName', var_name);

    xlim(ax, [0, 48]);
    ylim(ax, [0, max(var_48h) * 1.1]);
    xticks(ax, [0, 12, 24, 36, 48]);
    xticklabels(ax, {'7am', '7pm', '7am', '7pm', '7am'});
    xlabel(ax, 'Time (h)', 'FontSize', 12);
    ylabel(ax, var_name, 'FontSize', 12);
    title(ax, [var_name, ' over 48 hours'], 'FontSize', 14);
    grid(ax, 'on'); ax.GridAlpha = 0.3;
    legend(ax, 'show', 'FontSize', 11);
end
