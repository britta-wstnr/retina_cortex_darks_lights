% Plot ITC peaks (Fig. 6)

project_settings;

%% peaks plotting, identified with plot_itc_traces.m

on_retina = [
    NaN, 0.022, 0.024, NaN, 0.029, 0.034, 0.020, 0.039, 0.050, 0.066;
    0.048, 0.028, 0.034, 0.026, NaN, 0.032, 0.027, 0.037, 0.035, 0.048;
    0.029 0.027 NaN 0.027 0.025 0.026 0.037 0.026 0.024 0.031;
    0.013 0.029 0.029 0.022 0.030 0.027 0.023 0.027 0.027 0.025
    ] * 1000;

on_brain = [
    0.081, 0.069, 0.053, 0.058, 0.071, 0.074, 0.056, 0.047, NaN, 0.114;
    0.070, 0.066, 0.059, 0.096, 0.065, 0.084, 0.096, 0.046, 0.073, 0.087;
    0.055, 0.078, 0.065, 0.044, 0.084, 0.079, 0.082, 0.053, 0.076, 0.089;
    0.045, 0.071, 0.058, 0.057, 0.059, 0.054, 0.068, 0.051, NaN, 0.071;
    ] * 1000;

off_retina = [
    0.034, NaN, NaN, 0.040, 0.029, 0.033, 0.045, 0.026, NaN, 0.039] * 1000;

off_brain = [
    0.055, 0.059, 0.052, 0.047, 0.045, 0.070, 0.063, 0.067, 0.043, 0.065
            ] * 1000;


%% plot the median

colors = [
    239, 138,  98
    178,  24,  43
    103, 169, 207
     33, 102, 172
         ] / 255;
grey = [.7, .7, .7];
light_grey = [.9, .9, .9];

my_f = figure;
hold all;
y_axes = [1, 2.16, 3, 4];
x_min = 0; x_max = 100;

% plot patches
p_width = 0.3;
for ii = 1:length(y_axes)
    patch([x_min, x_max, x_max, x_min], ...
          [ii - p_width, ii - p_width, ii + p_width, ii + p_width], ...
          light_grey, 'EdgeColor', light_grey);
end


% plot the averages
p3 = errorbar(nanmedian(on_brain, 2), y_axes, iqr(on_brain, 2) ./ 2, ...
              'horizontal', ...
              'markeredgecolor', colors(3, :), 'markerfacecolor', ...
              colors(3, :), 'marker', 's', 'linestyle', 'none', ...
              'markersize', 14); hold on
p4 = errorbar(nanmedian(off_brain, 2), 1.84, iqr(off_brain) ./ 2, ...
              'horizontal', ...
              'markeredgecolor', colors(4, :), ...
              'markerfacecolor', colors(4, :), 'marker', 'o', 'linestyle', ...
              'none', 'markersize', 11); hold on
p1 = errorbar(nanmedian(on_retina, 2), y_axes, iqr(on_retina, 2) ./ 2, ...
              'horizontal', 'markeredgecolor', colors(1, :), ...
              'markerfacecolor', colors(1, :), 'marker', 's', 'linestyle', ...
              'none', 'markersize', 14); hold on
p2 = errorbar(nanmedian(off_retina, 2), 1.84, iqr(off_retina) ./ 2, ...
              'horizontal', 'markeredgecolor', colors(2, :), ...
              'markerfacecolor', colors(2, :), 'marker', 'o', 'linestyle', ...
              'none', 'markersize', 11); hold on

p1.Color = [0.5, 0.5, 0.5];
p1.CapSize = 5;
p1.LineWidth = 1.5;

p2.Color = [0.5, 0.5, 0.5];
p2.CapSize = 5;
p2.LineWidth = 1.5;

p3.Color = [0.5, 0.5, 0.5];
p3.CapSize = 5;
p3.LineWidth = 1.5;

p4.Color = [0.5, 0.5, 0.5];
p4.CapSize = 5;
p4.LineWidth = 1.5;

% x axis
set(gca, 'xlim', [x_min, x_max]);
xlabel('Time (ms)')

% y axis
set(gca, 'ylim', [0, 6]);
ax = gca;
ax.YTick = [1:4];
ax.YTickLabel = {'55-75', '75-95', '105-125', '125-145'};
ylabel('Frequency (Hz)');
ax.YGrid = 'on';

% legend
legend([p1, p2, p3, p4], {'retina, lights ON', 'retina, lights OFF', ...
        'visual cortex, lights ON', 'visual cortex, lights OFF'}, ...
       'location', 'northwest');
legend('boxoff');

% other figure characeristics
set(gca, 'fontsize', 18);
set(gcf, 'color', [1, 1, 1]);
set(my_f, 'position', [100, 100, 1049, 700]);

print(my_f, '-dpdf', '-bestfit', '-r1000', fullfile(fig_dir, ...
                                                  'peak_times.pdf'));
