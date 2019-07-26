% Plots the retinal evoked potentials (Fig. 1)

project_settings;

%% Load the data

for nn = 1:length(subjs)
    in_dir = fullfile(base_dir, subjs{nn});

    load(fullfile(in_dir, data_on_erg_fname));
    dataon{nn} = dataon_clean2;

    load(fullfile(in_dir, data_off_erg_fname));
    dataoff{nn} = dataoff_clean2;

end

%% Find peak times, ON data
baseline_window = [-.15 -.01];

for nn = 1:length(subjs)

    cfg = [];
    cfg.baseline_window = [-.15 -0.05];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    dataon_detrend{nn} = ft_preprocessing(cfg, dataon{nn})

    tlk_on{nn} = ft_timelockanalysis([], dataon_detrend{nn});

    if nn == 1
        atb = dsearchn(tlk_on{nn}.time', 0);
        ztb = dsearchn(tlk_on{nn}.time', 0.25);

        alltime = tlk_on{nn}.time(atb:ztb);
    end

    [~, retp] = findpeaks(squeeze(tlk_on{nn}.avg(2, atb:ztb)), ...
                          'minpeakwidth', 5, 'npeaks', 3);
    [~, retpmin] = findpeaks(squeeze(tlk_on{nn}.avg(2, atb:ztb)*-1), ...
                             'minpeakwidth', 5, 'npeaks', 2);

    rett = alltime(retp);
    rette = sprintf('%0.3f  ', rett);
    rettmin = alltime(retpmin);
    rettemin = sprintf('%0.3f  ', rettmin);

    figure
    plot(tlk_on{nn}.time, tlk_on{nn}.avg(2, :), ...
         'color', colors_rb(1, :), 'linewidth', 3);

    title(sprintf('%s, peak: %s, trough: %s', subjs{nn}, rette, rettemin))
    set(gca, 'xlim', [-.25 .25])
end

%% Find peak times, OFF data

for nn = 1:length(subjs)

    cfg = [];
    cfg.baseline_window = [-0.15 -0.05];
    cfg.detrend = 'yes'
    cfg.demean = 'yes'
    dataoff_detrend{nn} = ft_preprocessing(cfg, dataoff{nn})

    tlk_off{nn} = ft_timelockanalysis([], dataoff_detrend{nn});

    if nn == 1
        atb = dsearchn(tlk_off{nn}.time', 0);
        ztb = dsearchn(tlk_off{nn}.time', 0.1);

        alltime = tlk_on{nn}.time(atb:ztb);
    end

    [~,retp] = findpeaks(squeeze(tlk_on{nn}.avg(2, atb:ztb)), ...
                         'npeaks', 3);

    rett = alltime(retp);
    rette = sprintf('%0.3f  ', rett);

    figure
    plot(tlk_off{nn}.time, tlk_off{nn}.avg(2, :), ...
         'color', colors_rb(1, :), 'linewidth', 3);

    title(sprintf('%s, peak: %s', subjs{nn}, rette))
    set(gca, 'xlim', [-.15 .15])
end

%% Plot averaged light onset evoked potential

tlk_mean_on = ft_timelockgrandaverage([], tlk_on{:})

h = figure;
plot(tlk_mean_on.time * 1000, [tlk_mean_on.avg(2, :)] * 1e6, ...
     'color', colors_rb(3, :), 'linewidth', 3);
set(gca, 'xlim', [-100 300])
ylims = ylim
line([0 0], [ylims], 'color', [0 0 0], 'linewidth', 1.5);
set(gca, 'ylim', [ylims])

set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)

curr_ax = gca;
fix_plot(curr_ax, h);

ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

set(h, 'Position', [100, 100, 1049, 895]);
print(h,'-dpdf', '-bestfit', fullfile(fig_dir, 'ev_on.pdf'))

%% Plot averaged light offset evoked potential

tlk_mean_off = ft_timelockgrandaverage([], tlk_off{:})

h = figure;
plot(tlk_mean_off.time * 1000, [tlk_mean_off.avg(2, :)] * 1e6, ...
     'color', colors_rb(4, :), 'linewidth', 3);
set(gca, 'xlim', [-100 100])
ylims = ylim
line([0 0],[ylims], 'color', [0 0 0], 'linewidth', 1.5);
set(gca, 'ylim', [ylims])

set(gcf, 'color', [1 1 1])
set(gca,'FontSize',18)

curr_ax = gca;
fix_plot(curr_ax, h);

ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

set(h, 'Position', [100, 100, 1049, 895]);
print(h,'-dpdf', '-bestfit', fullfile(fig_dir, 'ev_off.pdf'))
