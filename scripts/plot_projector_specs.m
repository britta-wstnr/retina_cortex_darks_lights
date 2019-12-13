% Plot projector raise and fall times

project_settings;

%% Load raw data file

nn = 1;
proc_dir = 'long_subj_dir_identifier'

meg_file_base = ['proj_identifier', subjs{nn}];
meg_file_ending = '.fif';
meg_file_num = {'001', '002', '003'};

time = [];
trial = [];

% concatenate MEG files
for jj = 1:3

    meg_file = fullfile(proc_dir, ...
                        [meg_file_base, meg_file_num{jj}, ...
                        meg_file_ending]);

    cfg = [];
    cfg.dataset = meg_file;
    cfg.demean = 'yes';
    datafull = ft_preprocessing(cfg);
    time = [time, datafull.time{1}];
    trial = [trial, datafull.trial{1}];
end

datafull.time{1}  = time;
datafull.trial{1} = trial;

%% Select photodiode and define trials

clear cfg
cfg.channel = {'MISC001'};
photodiode = ft_selectdata(cfg, datafull);

photodiode.trial{1} = smooth(photodiode.trial{1}, 20)';

% ON data
indices = photodiode.trial{1} >= 0.02;
strindices = num2str(indices,'%-d');
onind = strfind(strindices, '01');
onind = onind + 1;
% correct onind for presentation screen at the end
onind = onind(1:250);

% OFF data
offind = strfind(strindices, '10');


%% Cut data in trials

clear cfgtrig
cfgtrig.trl = [onind' - 1000 ...
               onind' + 1000 ...
               repmat(-1000, ...
                      length(onind),1) repmat(1, length(onind),1)];
data_on = ft_redefinetrial(cfgtrig, datafull);

cfg = [];
cfg.channel = {'MISC001'};  % we only need the photo diode here
data_on = ft_selectdata(cfg, data_on);

clear cfgtrig
cfgtrig.trl = [offind' - 1000 ...
               offind' + 1000 ...
               repmat(-1000, ...
                      length(offind), 1) repmat(1, length(offind), 1)];
data_off = ft_redefinetrial(cfgtrig, datafull);

cfg = [];
cfg.channel = {'MISC001'};  % we only need the photo diode here
data_off = ft_selectdata(cfg, data_off);


%% Plot photodiode rise and fall times

g = figure;
hold all
start = dsearchn(data_on.time{1}', -0.01);
stop = dsearchn(data_on.time{1}', 0.01);
for ii = 1:100
    plot(data_on.time{ii}(start:stop).* 1000, ...
        data_on.trial{ii}(:, start:stop), '-', ...
        'color', colors_rb(4,:), ...
        'markeredgecolor', colors_rb(3, :), 'markerfacecolor', colors_rb(3, :));
end
set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)
ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

% properties of figure
curr_ax = gca;
fix_plot(curr_ax, g);

% save
print(g,'-dpdf', '-bestfit',fullfile(fig_dir, 'rise_times.pdf'))


h = figure;
hold all
start = dsearchn(data_off.time{1}', -0.01);
stop = dsearchn(data_off.time{1}', 0.01);
for ii = 1:100
    plot(data_off.time{ii}(start:stop) .* 1000, ...
        data_off.trial{ii}(:, start:stop), '-', ...
        'color', colors_rb(4, :), ...
        'markeredgecolor', colors_rb(3, :), 'markerfacecolor', colors_rb(3, :));
end
set(gcf, 'color', [1 1 1])
set(gca, 'FontSize', 18)
ylabel('Amplitude (\muV)')
xlabel('Time (ms)')

% properties of figure
curr_ax = gca;
fix_plot(curr_ax, h);

% save
print(h,'-dpdf', '-bestfit',fullfile(fig_dir, 'fall_times.pdf'))
