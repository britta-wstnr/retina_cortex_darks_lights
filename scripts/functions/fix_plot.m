function fix_plot(current_axes, figure)
   % Prettify MATLAB plots

    % set box property to off and remove background color
    set(current_axes, 'box', 'off', 'color', 'none')
    b = axes('Position', get(current_axes, 'Position'), 'box', 'on', ...
             'xtick', [], 'ytick', []);  % create new, empty axes with box
                                         % but without ticks
    axes(current_axes)  % set original axes as active
    linkaxes([current_axes b])  % link axes in case of zooming

    % adjust size
    set(figure, 'Position', [100, 100, 1049, 895]);