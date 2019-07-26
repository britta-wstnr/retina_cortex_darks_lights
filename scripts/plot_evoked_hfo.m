% Plot evoked retinal high frequency oscillations (Fig. 2)

project_settings;

%% load the data
% this results in one cell with 5 data structures (freq-bands), within each
% data structure are 10 trials (= subjects)
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
for nn = 1:length(subjs)
    tlk_on_erg{nn} = filter_osc_pot(data_on{nn}, 'on', 'yes', filter_specs_on);
    tlk_off_erg{nn} = filter_osc_pot(data_off{nn}, 'off', 'yes', ...
                                     filter_specs_off);
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


%% Plot oscillatory potential for light ON

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
clear data_on, data_off;

%% load the data - MEG

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

for stim = 1:2
    stim_def = stim_defs{stim};
    for nn = 1:length(subjs)
        if strcmp(stim_def, 'on')
            data_stim = filter_osc_pot(data_on{nn}, 'on', 'no', ...
                                       filter_specs_on);
        elseif strcmp(stim_def, 'off')
            data_stim = filter_osc_pot(data_off{nn}, 'off', 'no', ...
                                       filter_specs_off);
        end

        cfg = [];
        cfg.latency = [-0.15, 0.25];
        data_stim = ft_selectdata(cfg, data_stim);
        % run the spatial filter on all the data
        run_beamformer;

        % get the maximum voxel based on active time window:
        t_int = dsearchn(source_out.time', [0, 0.15]');
        pow = source_out.avg.pow;
        for ii = 1:length(source_out.inside)
            if ~isempty(source_out.avg.mom{ii})
                pow(ii) = mean(source_out.avg.mom{ii}(t_int(1):t_int(2)).^2);
            end
        end
        [val, idx] = max(pow);
        max_pos(stim, nn, :) = source_out.pos(idx, :);


        if strcmp(stim_def, 'on')
            tlk_on_meg{nn}.time = source_out.time;
            tlk_on_meg{nn}.avg = source_out.avg.mom{idx};
        else
            tlk_off_meg{nn}.time = source_out.time;
            tlk_off_meg{nn}.avg = source_out.avg.mom{idx};
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
% t_shift = 0;  % to also have the original one


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
