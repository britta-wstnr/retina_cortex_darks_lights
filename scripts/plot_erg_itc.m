% Plots the ITC data for retinal light offset activity (Fig. 2)

project_settings;

baselinewindow = [-.15 -.05];

%% Read data

% this results in one cell with 5 data structures (freq bands), within each
% data structure are 10 trials (subjects)
for nn = 1:length(subjs)

    indir       = fullfile(basedir, subjs{nn}, 'hilbert_out');

    load(fullfile(indir, 'erg_OFF_rightstats.mat'));

    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.latency = [-.15 .15];
    cfg.channel = 'EOG002';
    data_tf = ft_selectdata(cfg, data_tf);

    for ff = 1:size(freq_bands, 1);
        data_subjs_on{ff}.time = cell(1, length(subjs));
        data_subjs_on{ff}.time(:) = {data_tf.time};
        data_subjs_on{ff}.trialinfo = repmat(1,1,length(subjs));
        data_subjs_on{ff}.label = data_tf.label;

        data_subjs_on{ff}.trial{nn} = squeeze(data_tf.itc(:,ff,:))';
    end
end

data = data_subjs_on;

%% run group-level stats

for ii = 1:size(freq_bands, 1)
    nemo_subjectstats_sensorlevel
end

%% assembles data_tested{:} into a composite source_tf structure
data_plot = data_tested{1};
data_plot.avg = zeros(numel(data_tested{1}.label), size(freq_bands,1), ...
                      size(data_tested{1}.time, 2));
data_plot.stat = data_plot.avg;
data_plot.pval = data_plot.avg;

data_plot.freq_bands = freq_bands;

% This is needed if output should be plotted with fieldtrip's plotting
% functions - ft_singleplotTFR or ft_multiplotTFR
data_plot.freq = mean(data_plot.freq_bands');

for jj = 1:size(freq_bands, 1)
    data_plot.avg(:, jj, :) = reshape(data_tested{jj}.avg, ..
                                      numel(data_tested{1}.label), 1, ...
                                      size(data_tested{1}.time, 2));
    % data_plot.itc(:, jj, :) = reshape(data_tested{jj}.itc, ...
    %                                 numel(data_tested{1}.label), 1, ...
    %                                 size(data_tested{1}.time, 2));

    data_plot.stat(:, jj, :) = reshape(data_tested{jj}.stat, ...
                                       numel(data_tested{1}.label), 1, ...
                                       size(data_tested{1}.time, 2));
    data_plot.pval(:, jj, :) = reshape(data_tested{jj}.pval, ...
                                       numel(data_tested{1}.label), 1, ...
                                       size(data_tested{1}.time, 2));

end

data_plot.dimord = 'chan_freq_time';

%% enable different MCP corrections:

testcorrect = 'fdr'

switch testcorrect
    case 'fdr'
        addpath /fdr_bh/

        pvals = 10 .^ (data_plot.pval);
        fdr = fdr_bh(pvals, 0.05, 'pdep', 'yes');
        data_plot.mask = zeros(size(data_plot.pval));
        data_plot.mask = logical(fdr);
    case 'fdr_yek'
        addpath /fdr_bh/

        pvals = 10 .^ (data_plot.pval);
        fdr = fdr_bh(pvals, 0.05, 'dep', 'yes');
        data_plot.mask = zeros(size(data_plot.pval));
        data_plot.mask = logical(fdr);

    case 'bonf'
        % bonferroni correction and conversion to log scale
        plim=log10(0.05/prod(size(data_plot.pval)));
        mask = zeros(size(data_plot.pval));
        idx = find(data_plot.pval < plim);
        mask(idx) = 1;
        data_plot.mask = logical(mask);
end


%% plot
h=figure;

cfg = [];
cfg.funparameter = 'stat';
cfg.maskparameter = 'mask';
cfg.marking = 1;
cfg.markingclust = 5;
cfg.zeroline = 1;
cfg.channel = 'EOG002';

nemo_plot_hilbertchannel(cfg, data_plot);

set(gca, 'fontsize', 18)
title(' ')
set(gcf, 'color', [1 1 1])

%% save figure

fig_dir = fullfile(base_dir, 'figures_paper')

print(h,'-dpdf', '-bestfit', '-opengl',...
      fullfile(fig_dir, 'itc_erg_off_broadb.pdf'))
