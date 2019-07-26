% Plot ITC timecourses (Fig. 5) and identifying ITC peak times


project_settings;
hilbert_dir = fullfile(base_dir, 'hilbert_results');

%% get occipital voxels

load(fullfile(base_dir, 'atlasgrid.mat'));   % same positions as source_tf
search_grid = {'Occipital_Mid_L', 'Occipital_Sup_L', 'Occipital_Inf_L', ...
              'Occipital_Mid_R', 'Occipital_Inf_R', 'Occipital_Sup_R', ...
              'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R'};

atlas_reg= [];
for ii = 1:length(search_grid)
    dummy = strmatch(search_grid{ii}, atlasgrid.tissuelabel);
    dummy = find(atlasgrid.gridlabel == dummy);
    atlas_reg = [atlas_reg, dummy'];
end
atlas_reg = unique(atlas_reg);  % already does the sorting


%% Process light ON or OFF data

% which to plot, light on or off:
plot_it = 'OFF';

% load subject data
% this results in one cell with 5 data structures (freq bands), within each
% data structure are 10 trials (subjects)
for nn = 1:length(subjs)

    % indivudal directory
    in_dir = fullfile(base_dir, subjs{nn}, 'hilbert_out');

    % RETINA data
    load(fullfile(in_dir, ['erg_', plot_it, '_rightstats.mat']));

    rt_min = dsearchn(data_tf.time', -.15);
    rt_max = dsearchn(data_tf.time', .15);

    for ii = 1:4  % four frequency bands
        retina_itc(ii,nn,:)  = data_tf.itc(2, ii, rt_min:rt_max);
    end

    % BRAIN data
    load(fullfile(in_dir, ['sourcetf_medcov_', plot_it, '.mat']));
    inside_idx = find(source_tf.inside);

    bt_min = dsearchn(source_tf.time', -.15);
    bt_max = dsearchn(source_tf.time', .15);
    peak_tb = dsearchn(source_tf.time', .02);

    for ii=1:4  % four relevant frequency bands
        for vox = 1:length(atlas_reg)
            % find maximum 0.02 - 0.15 sec
            [maxval(vox)] = max(source_tf.avg.itc{atlas_reg(vox)}(ii, ...
                                                              peak_tb:bt_max));
        end

        [~, maxvox(ii, nn)] = max(maxval);
        brain_itc(ii, nn, :) = source_tf.avg.itc{atlas_reg(maxvox(ii, ...
                                                          nn))}(ii, ...
                                                          bt_min:bt_max);
    end
end

%% Identify peak times

% preallocation
max_ret = []; max_ret_idx = []; max_ret_tp= [];
max_brain = []; max_brain_idx = []; max_brain_tp = [];

% time points
tp_a = dsearchn(source_tf.time', -.15);
tp_b = dsearchn(source_tf.time', .15);
tp_c = dsearchn(source_tf.time(tp_a:tp_b)', 0);
tp_d = dsearchn(source_tf.time(tp_a:tp_b)', 0.1);

my_time = source_tf.time(tp_a:tp_b);
my_time_2 = source_tf.time(tp_c:tp_d);

for ii = 1:length(subjs)
    h = figure
    for ff = 1:3  % freqbands, change to needs
        subplot(1, 2, ff-1)
        plot(source_tf.time(tp_a:tp_b), squeeze(retina_itc(ff, ii, :)), ...
             'color', colors_rb(2, :), 'linewidth', 3); hold on;
        plot(source_tf.time(tp_a:tp_b), squeeze(brain_itc(ff, ii, :)), ...
             'color', colors_rb(4,:), 'linewidth', 3); hold on;

        % identify the peaks
        [~,ret_p] = findpeaks(squeeze(retina_itc(ff, ii, :)), ...
                             'minpeakwidth', 5);
        [~, brain_p] = findpeaks(squeeze(brain_itc(ff, ii, :)), ...
                                'minpeakwidth', 5);

        ret_t = my_time(ret_p);
        ret_t = ret_t(find(ret_t > 0));
        ret_t = ret_t(find(ret_t < .1));

        brain_t = my_time(brain_p);
        brain_t = brain_t(find(brain_t > 0));
        brain_t = brain_t(find(brain_t < .1));

        set(gca, 'xlim', [-Inf Inf])

        ret_t_string = sprintf('%0.3f  ', ret_t);
        brain_t_string = sprintf('%0.3f  ', brain_t);
        title(sprintf('%d - %d Hz; peaks retina %s; peaks brain: %s', ...
                      freqbands(ff, 1), freqbands(ff, 2), ...
                      ret_t_string, brain_t_string));

        [max_ret(ii, ff), max_ret_idx(ii, ff)] = ...
            max(squeeze(retina_itc(ff, ii, tp_c:tp_d)));
        [max_brain(ii, ff), max_brain_idx(ii, ff)] = ...
            max(squeeze(brain_itc(ff, ii, tp_c:tp_d)));
        max_ret_tp(ii, ff) = my_time_2(max_ret_idx(ii, ff));
        max_brain_tp(ii, ff) = my_time_2(max_brain_idx(ii, ff));

        if(ff==2)
            legend({'retina', 'brain'}, 'location', 'northwest')
        end

    end
    suptitle(sprintf('ITC %s, subject %s', plot_it, subjs{ii}))
    set(gcf, 'color', [1,1,1]);
    set(h, 'Position', [100, 100, 1049, 895]);

end


%% plot ITC subjects (median and dispersion)

pt_min = dsearchn(source_tf.time(bt_min:bt_max)', 0);
pt_max = dsearchn(source_tf.time(bt_min:bt_max)', 0.15);

all_time = source_tf.time(bt_min:bt_max);
peak_search_time = source_tf.time(pt_min:pt_max);

if strcmp(plot_it, 'ON')
    color_brain = colors_rb(1, :);
    color_retina = colors_rb(3, :);
else
    color_brain = colors_rb(2, :);
    color_retina = colors_rb(4, :);
end


for ff = 1:4  % freq_bands
    h = figure

    [l(1), p] = boundedline(source_tf.time(bt_min:bt_max), ...
                            median(squeeze(retina_itc(ff, :, :))), ...
                            iqr(squeeze(retina_itc(ff, :, :))) / 2, ...
                            'alpha', 'g', 'transparency', 0.2, ...
                            'cmap', color_retina);
    hold on;
    [l(2), p] = boundedline(source_tf.time(bt_min:bt_max), ...
                            median(squeeze(brain_itc(ff, :, :))), ...
                            iqr(squeeze(brain_itc(ff, :, :))) / 2,  ...
                            'alpha', 'g', 'transparency', 0.2, ...
                            'cmap', color_brain);
    hold on;
    set(l(1), 'LineWidth', 3);
    set(l(2), 'LineWidth', 3);
    title(sprintf('%d - %d Hz', freq_bands(ff, 1), freq_bands(ff, 2)));
    set(gca, 'ylim', [0, 1.0])
    set(gca, 'YTick', [0:0.2:1.0]);
    set(gca, 'FontSize', 18)
    set(gca, 'xlim', [-Inf Inf])

    if(ff == 1)
        legend([l], {'Retina', 'Occipital cortex'}, 'location', 'northwest')
    end

    set(gcf, 'color', [1, 1, 1])
    set(h, 'Position', [100, 100, 1049, 895]);

    print(h,'-dpdf', '-bestfit', '-r1000', ...
          fullfile(fig_dir, ['itc_', plot_it, '_',  ...
                        num2str([freq_bands(ff)]), '.pdf']))
end
