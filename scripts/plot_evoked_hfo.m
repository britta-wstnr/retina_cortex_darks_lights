% Plot evoked retinal high frequency oscillations (Fig. 2)
% Compute and plot Granger causality results (Fig. 6)

project_settings;

%% load the data
% this results in one cell with 5 data structures (freq-bands), within each
% data structure are 10 trials (= subjects)

data_on = cell(10, 1);
data_off = cell(10, 1);

for nn = 1:length(subjs)
    in_dir = fullfile(base_dir, subjs{nn});
    load(fullfile(in_dir, data_on_erg_fname));
    data_on{nn} = dataon_clean2;

    load(fullfile(in_dir, data_off_erg_fname));
    data_off{nn} = dataoff_clean2;

end

%% process the data - ERG
% filtering the data and averaging across trials per subject
% light ON:

% preallocation:
tlk_on_erg = cell(10, 1);
on_erg_granger = cell(10, 1);
tlk_off_erg = cell(10, 1);
off_erg_granger = cell(10, 1);

for nn = 1:length(subjs)
    tlk_on_erg{nn} = filter_osc_pot(data_on{nn}, 'on', 'yes', filter_specs_on);
    on_erg_granger{nn} = filter_osc_pot(data_on{nn}, 'on', 'no', ...
        filter_specs_on);
    tlk_off_erg{nn} = filter_osc_pot(data_off{nn}, 'off', 'yes', ...
        filter_specs_off);
    % only high-pass filter for later granger estimation
    off_erg_granger{nn} = filter_osc_pot(data_off{nn}, 'on', 'no', ...
        filter_specs_on);
end

%% Plot oscillatory potential for light ON

% average across subjects:
tlk_mean_on = ft_timelockgrandaverage([], tlk_on_erg{:});

% plot and save figure
h = figure;
plot(tlk_mean_on.time * 1000, ...  % milliseconds
    [tlk_mean_on.avg(2, :)] * 1e6, ...
    'color', colors_rb(3, :), 'linewidth', 3);
line([0 0], [-5 5], 'color', [0 0 0], 'linewidth', 1.5);  % line at zero
set(gca, 'ylim', [-5 5])
set(gca, 'xlim', [-150 250])
set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)
ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

% properties of figure
curr_ax = gca;
fix_plot(curr_ax, h);

% save
print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'ev_hfo_on.pdf'))


%% Plot oscillatory potential for light OFF

% average across subjects:
tlk_mean_off = ft_timelockgrandaverage([], tlk_off_erg{:});

% plot and save figure
h = figure;
plot(tlk_mean_off.time * 1000, ...  % milliseconds
    [tlk_mean_off.avg(2, :)] * 1e6, ...
    'color', colors_rb(3, :), 'linewidth', 3);
line([0 0], [-1 1], 'color', [0 0 0], 'linewidth', 1.5);  % line at zero
set(gca, 'ylim', [-1 1])
set(gca, 'xlim', [-150 250])
set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)
ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

% properties of figure
curr_ax = gca;
fix_plot(curr_ax, h);

% save
print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'ev_hfo_off.pdf'))

%% memory management
clear data_on data_off;

%% load the data - MEG

% preallocation:
data_on = cell(10, 1);
data_off = cell(10, 1);
for nn = 1:length(subjs)
    in_dir = fullfile(base_dir, subjs{nn});
    load(fullfile(in_dir, data_on_meg_fname));
    load(fullfile(in_dir, data_off_meg_fname));

    % cut to the needed channels to save time
    cfg = [];
    cfg.channel = meg_channel_sel;
    dataon_clean = ft_selectdata(cfg, dataon_clean);
    dataoff_clean = ft_selectdata(cfg, dataoff_clean);

    data_on{nn} = dataon_clean;
    data_off{nn} = dataoff_clean;

end

%% process the data - MEG
% filtering the data, source reconstructing, and averaging across trials
% per subject

% preallocation
tlk_on_meg = cell(10, 1);
tlk_on_meg_granger = cell(10, 1);
tlk_off_meg = cell(10, 1);
off_meg_granger = cell(10, 1);
on_meg_granger = cell(10, 1);
max_pos = zeros(2, 10, 3);
indices = zeros(2, 10);

for stim = 1:4  % also account for special granger data

    for nn = 1:length(subjs)
        if stim == 1 || stim == 4
            data_stim = filter_osc_pot(data_on{nn}, 'on', 'no', ...
                filter_specs_on);
        elseif stim == 2
            data_stim = filter_osc_pot(data_off{nn}, 'off', 'no', ...
                filter_specs_off);
        elseif stim == 3
            % get high-pass filtered data for later granger
            data_stim = filter_osc_pot(data_off{nn}, 'on', 'no', ...
                filter_specs_on);
        end

        cfg = [];
        cfg.latency = [-0.15, 0.25];
        data_stim = ft_selectdata(cfg, data_stim);
        % run the spatial filter on all the data
        run_beamformer;

        % get the maximum voxel based on active time window:
        if stim == 1 || 2
            t_int = dsearchn(source_out.time', [0, 0.15]');
            pow = source_out.avg.pow;
            for ii = 1:length(source_out.inside)
                if ~isempty(source_out.avg.mom{ii})
                    pow(ii) = mean(source_out.avg.mom{ii}(t_int(1):t_int(2)).^2);
                end
            end
            [val, idx] = max(pow);
            indices(stim, nn) = idx;
            max_pos(stim, nn, :) = source_out.pos(idx, :);
        end


        if stim == 1
            tlk_on_meg{nn}.time = source_out.time;
            tlk_on_meg{nn}.avg = source_out.avg.mom{idx};
            % keep some extra info around
            tlk_on_meg{nn}.fsample = data_on{1}.fsample;
            tlk_on_meg{nn}.label = 'max_vox';
        elseif stim == 2
            tlk_off_meg{nn}.time = source_out.time;
            tlk_off_meg{nn}.avg = source_out.avg.mom{idx};
            % keep some extra info around
            tlk_off_meg{nn}.fsample = data_on{1}.fsample;
            tlk_off_meg{nn}.label = 'max_vox';
        elseif stim == 3 || stim == 4
            if stim == 3
                cond_equals = 2;  % treat like OFF data
            else
                cond_equals = 1;  % treat like ON data
            end

            % we need all trials here, not just the average
            spat_filt = source_out.avg.filter{indices(cond_equals, nn)};

            out_tmp = cell(length(data_stim.trial), 1);
            out_time = out_tmp;
            for tr = 1:length(data_stim.trial)
                out_tmp{tr} = spat_filt * data_stim.trial{tr};
                out_time{tr} = source_out.time;
            end

            if stim == 3
                % store in output structure:
                off_meg_granger{nn}.time = out_time;
                off_meg_granger{nn}.trial = out_tmp;
                % keep some extra info around
                off_meg_granger{nn}.fsample = data_on{1}.fsample;
                off_meg_granger{nn}.label = 'max_vox';
            else
                on_meg_granger{nn}.time = out_time;
                on_meg_granger{nn}.trial = out_tmp;
                % keep some extra info around
                on_meg_granger{nn}.fsample = data_on{1}.fsample;
                on_meg_granger{nn}.label = 'max_vox';
            end
        end
    end
end


%% Plot oscillatory potential and cortical high frequency activity

% indicate subject here; cleanest one is 10
nn = 10;

t_int = dsearchn(tlk_on_erg{nn}.time', [-0.15, 0.25]');
[~, t_ret] = max(abs(tlk_on_erg{nn}.avg(2, t_int(1):t_int(2))));
[~, t_cor] = max(abs(tlk_on_meg{nn}.avg));

% calculate shift in time
shift = t_cor - t_ret;
t_shift = shift / 1500;  % sampling frequency
t_shift = 0;  % to also have the original one


% plot and save figure
h = figure;
plot(tlk_on_erg{nn}.time * 1000, ...  % milliseconds
    [tlk_on_erg{nn}.avg(2, :)] * 1e6, ...
    'color', colors_rb(3, :), 'linewidth', 3); hold on
plot((tlk_on_meg{nn}.time - t_shift) * 1000, ...  % milliseconds
    [tlk_on_meg{nn}.avg] *-1 , ...
    'color', colors_rb(1, :), 'linewidth', 3);

line([0 0], [-5.5 5.5], 'color', [0 0 0], 'linewidth', 1.5);  % line at zero
set(gca, 'ylim', [-5.5 5.5])
set(gca, 'xlim', [-150 250])
set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)
ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

% % properties of figure
% curr_ax = gca;
% fix_plot(curr_ax, h);

print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'ev_hfo_compare_timeshift.pdf'))


%% Compute spectrally-resolved Granger causality

% preallocation
on_granger_cort_ret = cell(10, 1);
on_granger_ret_cort = cell(10, 1);
off_granger_cort_ret = cell(10, 1);
off_granger_ret_cort = cell(10, 1);
on_flip_cort_ret = cell(10, 1);
on_flip_ret_cort = cell(10, 1);
off_flip_cort_ret = cell(10, 1);
off_flip_ret_cort = cell(10, 1);


for flips = 1:2
    flip_it = flips - 1;

    for stim = 1:2
        stim_def =  stim_defs{stim};

        if strcmp(stim_def, 'on')
            tlk_erg = on_erg_granger;
            tlk_meg = on_meg_granger;
        else
            tlk_erg = off_erg_granger;
            tlk_meg = off_meg_granger;
        end

        for nn = 1:length(subjs)
            % cut the data to the same length
            cfg = [];
            cfg.latency =  [0, 0.25];  %[-0.15, 0.25];
            granger_in = ft_selectdata(cfg, tlk_erg{nn});
            granger_in_meg = ft_selectdata(cfg, tlk_meg{nn});

            % substitue the first channel with the MEG data
            % the first channel is unused in the contrast anyway
            for tr = 1:length(granger_in.trial)
                granger_in.trial{tr}(1, :) = granger_in_meg.trial{tr};

                % reverse the spectra in the flipped cases
                if flip_it
                    granger_in.trial{tr}(1, :) = ...
                        fliplr(granger_in.trial{tr}(1, :));
                    granger_in.trial{tr}(2, :) = ...
                        fliplr(granger_in.trial{tr}(2, :));
                end
            end

            granger_in.label{1} = 'cortex';
            granger_in.label{2} = 'retina';

            % compute fourier analysis on data
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.output = 'fourier';
            cfg.foilim = [0 195];
            cfg.tapsmofrq = 5;
            cfg.pad = 3;
            freq = ft_freqanalysis(cfg, granger_in);

            % compute connectivity
            cfg = [];
            cfg.method = 'granger';
            granger_out = ft_connectivityanalysis(cfg, freq);

            % store as timelockdata for stats: both directions
            % make dummy timelock structure
            tmp.label = {'dummy'};
            tmp.dimord = 'chan_time';
            tmp.time = granger_out.freq;

            granger_ret_cort = tmp;  % retina to cortex
            granger_cort_ret = tmp;  % cortex to retina

            granger_cort_ret.avg(1, :) = granger_out.grangerspctrm(1, 2, :);
            granger_ret_cort.avg(1, :) = granger_out.grangerspctrm(2, 1, :);

            % store according to the stimulus:
            % make things very explicit to easily keep track
            if flip_it
                if strcmp(stim_def, 'on')
                    on_flip_cort_ret{nn} = granger_cort_ret;
                    on_flip_ret_cort{nn} = granger_ret_cort;
                else
                    off_flip_cort_ret{nn} = granger_cort_ret;
                    off_flip_ret_cort{nn} = granger_ret_cort;
                end
            else
                if strcmp(stim_def, 'on')
                    on_granger_cort_ret{nn} = granger_cort_ret;
                    on_granger_ret_cort{nn} = granger_ret_cort;
                else
                    off_granger_cort_ret{nn} = granger_cort_ret;
                    off_granger_ret_cort{nn} = granger_ret_cort;
                end
            end

        end  % loop over subjects
    end  % loop over stimulus conditions
end  % loop over flipped/non_flipped

close all;

%% Granger causality statistics

for stim = 1:2

    % handle the conditions
    if stim == 1
        granger_ret_cort = on_granger_ret_cort;
        granger_cort_ret = on_granger_cort_ret;
    else
        granger_ret_cort = off_granger_ret_cort;
        granger_cort_ret = off_granger_cort_ret;
    end

    % run the permutation test
    cfg = [];
    cfg.latency = [55 Inf];
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'
    cfg.alpha  = 0.025;
    cfg.correctm = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.clusterthreshold = 'nonparametric_individual';
    % cfg.correcttail = 'prob';
    cfg.numrandomization = 10000;

    n_sub = length(subjs);
    cfg.design(1, 1:2 * n_sub) = [ones(1, n_sub) 2 * ones(1, n_sub)];
    cfg.design(2, 1:2 * n_sub) = [1:n_sub 1:n_sub];
    cfg.ivar = 1; % the 1st row contains the independent variable
    cfg.uvar = 2; % the 2nd row contains the subject number

    stat{stim} = ft_timelockstatistics(cfg, granger_ret_cort{:}, ...
        granger_cort_ret{:})

end
%% Plot the statistical output and flipped time series

for flips = 1:2
    flip_it = flips - 1;

    for stim = 1:2
        stim_def = stim_defs{stim};

        % get the right time series for plotting
        if flip_it
            flip_str = 'reversed';
            flip_title = [', ', flip_str];  % title for plotting
            if stim == 1
                granger_ret_cort = on_flip_ret_cort;
                granger_cort_ret = on_flip_cort_ret;
            else
                granger_ret_cort = off_flip_ret_cort;
                granger_cort_ret = off_flip_cort_ret;
            end
        else
            flip_str = '';
            flip_title = '';
            if stim == 1
                granger_ret_cort = on_granger_ret_cort;
                granger_cort_ret = on_granger_cort_ret;
            else
                granger_ret_cort = off_granger_ret_cort;
                granger_cort_ret = off_granger_cort_ret;
            end
        end

        for n = 1:10
            plot_ret_cort(n, :) = granger_ret_cort{n}.avg;
            plot_cort_ret(n, :) = granger_cort_ret{n}.avg;
        end

        cut_off = dsearchn(granger_out.freq', 55);

        y_min = 0;
        y_max = 0.045;

        h = figure;

        % find the significant clusters if relevant
        if ~flip_it && sum(stat{stim}.mask) >= 1
            for ii = 1:size(stat{stim}.mask, 2) - 1;
                if ii == 1
                    if stat{stim}.mask(ii) == 1
                        cluster_start = ii;
                    end
                else
                    if stat{stim}.mask(ii-1) == 0 && stat{stim}.mask(ii) == 1
                        cluster_start = ii;
                    elseif stat{stim}.mask(ii) == 1 && ...
                            stat{stim}.mask(ii+1) == 0
                        cluster_stop = ii;
                    end
                end
            end

            % shade significant frequencies
            for ii = 1:length(cluster_start)
                x_min = stat{stim}.time(cluster_start(ii));
                x_max = stat{stim}.time(cluster_stop(ii));
                patch([x_min, x_max, x_max, x_min], ...
                    [y_min, y_min, y_max, y_max], ...
                    light_grey, 'EdgeColor', light_grey);
                display(stim_defs{stim})
                display(x_min)
                display(x_max)
            end
        end

        % plot the Granger causality estimates
        [l(1), p] = boundedline(granger_out.freq(cut_off:end), ...
            mean(plot_cort_ret(:, cut_off:end)), ...
            std(plot_cort_ret(:, cut_off:end)) / 2, ...
            'alpha', 'g', 'transparency', 0.2, ...
            'cmap', colors_granger(1, :));
        [l(2), p] = boundedline(granger_out.freq(cut_off:end), ...
            mean(plot_ret_cort(:, cut_off:end)), ...
            std(plot_ret_cort(:, cut_off:end)) / 2, ...
            'alpha', 'g', 'transparency', 0.2, ...
            'cmap', colors_granger(2, :));

        hold on;
        set(l(1), 'LineWidth', 3);
        set(l(2), 'LineWidth', 3);

        if strcmp(stim_defs{stim}, 'on')
            title(['Granger causality, light onset', flip_title]);
        else
            title(['Granger causality, light offset', flip_title]);
        end

        xlabel('Frequency (Hz)')
        ylabel('Granger causality');

        set(gca, 'ylim', [y_min, y_max])


        set(gca, 'YTick', [y_min:0.01:y_max]);
        %set(gca, 'FontSize', 18)
        set(gca, 'xlim', [-Inf Inf])

        legend([l], {'Cortex to retina', 'Retina to cortex'}, ...
            'location', 'northeast')


        set(gcf, 'color', [1, 1, 1])
        %set(h, 'Position', [100, 100, 1049, 895]);
        print(h,'-dpdf', '-bestfit',fullfile(fig_dir, ...
            ['granger_', stim_def, '_', flip_str,'.pdf']))

    end
end
