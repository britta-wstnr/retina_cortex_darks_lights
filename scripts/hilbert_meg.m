% Compute Hilbert transform of MEG data

% This script calls nemo_hilbertstats which uses parfor. Set parallel computing
% preferences before running.
% Also set stim_def to 'on' or 'off'

% script should be called by prep_hilbert_meg

project_settings;

%% specifications
ts_width = 6;
baseline_window = [-.15 -.05];
toi_lim  = [-.15 .25]        % this is needed after filtering

stim = stim_def ;
filter_type = 'firws'; % 'but' or 'firws'

% data already loaded. This is preprocessed + downsampled sensorlevel data.
if(strcmp(stim, 'on'))
    data = dataon_clean;
    display('Analyzing ON data');
elseif(strcmp(stim, 'off'))
    data = dataoff_clean;
    display('Analyzing OFF data');
else
    error('unknown stimulus condition!')
end

saving = 1;     % should data be saved?
plotting = 0;   % should results be plotted?

mediancov = 1;

% for plotting:
masking = 0;

%% run through filter bank, and obtain timelock and hilbert transform

% preallocation
timelock_bp = cell(1,size(freq_bands,1));
data_hilb  = cell(1,size(freq_bands,1));
filter_definitions = cell(1,size(freq_bands,1));

for ii = 1:size(freq_bands, 1)

    cfg=[];
    cfg.demean = 'yes';
    cfg.baseline_window = baseline_window;
    if(strcmp(filter_type, 'but'))
        % NB: default butterworth for quick testing; specify more advanced
        % filter for real analysis!
        cfg.bpfilter = 'yes';
        cfg.bpfreq = freq_bands(ii,:);
    elseif(strcmp(filter_type, 'firws'))
        cfg.hpfilter = 'yes';
        cfg.hpfreq = freq_bands(ii,1);
        cfg.hpfilttype = filter_type;
        cfg.hpfiltdf = ts_width;
        cfg.hpfiltdir = 'onepass-zerophase';
        cfg.hpfiltwintype = 'hann';
        cfg.plotfiltresp = 'no';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = freq_bands(ii,2);
        cfg.lpfilttype = filter_type;
        cfg.lpfiltdf = ts_width;
        cfg.lpfiltdir = 'onepass-zerophase';
        cfg.lpfiltwintype = 'hann';
    else
        error('I do not know the specified filter type, please add!')
    end

    data_bp = ft_preprocessing(cfg,data);

    % filtering first then snipping to shorter time interval avoids edge
    % artifacts
    cfg = [];
    cfg.toi_lim = toi_lim;
    data_bp = ft_redefinetrial(cfg, data_bp);

    cfgtl                  = [];
    cfgtl.covariance       = 'yes';
    cfgtl.covariancewindow = 'all';
    cfgtl.vartrllength     = 2;
    timelock_bp{ii}        = ft_timelockanalysis(cfgtl, data_bp);

    cfghilb = [];
    cfghilb.hilbert = 'complex';
    data_hilb{ii} = ft_preprocessing(cfghilb, data_bp);
end

% memory management
clear data_bp


%% median covariance computation

if(perfectori)
    baselinewindow_idx = dsearchn(data_hilb{1}.time{1}',baseline_window');
    activewindow_idx = dsearchn(data_hilb{1}.time{1}',activewindow');
end

if(mediancov)

    for ii = 1:size(freq_bands, 1)
        [med_cov, ~] = compute_median_cov(data_hilb{ii}, leadfield, ...
                                              meg_channel_sel)
        timelock_bp{ii}.cov = med_cov;
    end

%% create spatial filter using the non-Hilbert data

% prealloc
source_bp = cell(1,size(freq_bands,1));

for ii=1:size(freq_bands,1)
    cfg                   = [];
    cfg.method            = 'lcmv';
    cfg.grid              = leadfield;
    cfg.vol               = vol;
    cfg.keepfilter        = 'yes';
    cfg.lcmv.reducerank   = 'no';
    cfg.lcmv.fixedori     = 'yes'; % project on axis of most variance using SVD
    cfg.lcmv.projectnoise = 'yes';
    cfg.lcmv.weightnorm   = 'nai';
    cfg.lcmv.keepfilter   = 'yes';
    %    cfg.lcmv.lambda      = '10%';
    source_bp{ii}         = ft_sourceanalysis(cfg, timelock_bp{ii});
end


%% apply this spatial filter to complex-valued hilbert transform of the data

save_ram = 1;
tf_stats = 1;

disp(['Hilbert-sourcing started at ' datestr(now)])
t_start = tic;

for ii = 1:size(freq_bands,1)
    % change to ordinary "for" if parfor claims
    % "Attempt to serialize data which is too large."
    source_hilbtmp{ii} = ft_apply_spatialfilter(data_hilb{ii},source_bp{ii});
    data_hilb{ii} = [];

    % perform stats on Hilbert source trials
    nemo_hilbertstats
end

disp(['Hilbert-sourcing finished at ' datestr(now) ...
      ' and took ' num2str(toc(t_start)/60) ' minutes to run']);

if(save_ram)
    clear data_hilb
end

%% assembles source_hilb{:} into a composite source_tf structure

source_tf = source_hilb{1};
source_tf.freq_bands = freq_bands(1:end,:);
source_tf.freq = mean(source_tf.freq_bands');

for jj = 1:size(freq_bands, 1)
    for ii = 1:length(inside_idx)
        source_tf.avg.mom{inside_idx(ii)}(jj, :) = ...
            source_hilb{jj}.avg.mom{inside_idx(ii)};
        source_tf.avg.itc{inside_idx(ii)}(jj, :) = ...
            source_hilb{jj}.avg.itc{inside_idx(ii)};
        if(tf_stats)
            source_tf.stat{inside_idx(ii)}(jj, :) = ...
                source_hilb{jj}.stat{inside_idx(ii)};
            source_tf.pval{inside_idx(ii)}(jj, :) = ...
                source_hilb{jj}.pval{inside_idx(ii)};
        end
    end
end

source_tf.pos = template_grid.pos;  % supply MNI positions
source_tf.coordsys = 'mni';

% save data
if(strcmp(stim, 'on'))
    display('Saving ON data');
    save(fullfile(outdir, source_on_meg_fname), 'source_tf');
elseif(strcmp(stim, 'off'))
    display('Saving OFF data');
    save(fullfile(outdir, source_off_meg_fname), 'source_tf');
end

%% plotting

if(plotting)
    if(masking)
        %% plotting with thresholding (e.g. with p-value)

        cfg=[];
        cfg.mripath = t1_template_fname;
        cfg.funparameter = 'stat';
        cfg.maskparameter = 'msk';
        masktype = 'pval';
        maskthresh = log10(.01);

        if(isfield(source_tf,'dim'))
            source_tf = rmfield(source_tf,'dim');
        end

        source_tf.msk = zeros(size(source_tf.avg.mom,1), ...
                              length(source_tf.freq), ...
                              length(source_tf.time));

        switch(masktype)
            case 'itchigh'
                itc = reshape([source_tf.avg.itc{:}], ...
                              length(source_tf.freq), ...
                              length(source_tf.time), ...
                              length(inside_idx));
                itc = permute(itc, [3 1 2]);
                source_tf.msk(inside_idx, ...
                              1:length(source_tf.freq), :) = itc > maskthresh;
            case 'itclow'
                itc = reshape([source_tf.avg.itc{:}], ...
                              length(source_tf.freq), ...
                              length(source_tf.time), ...
                              length(inside_idx));
                itc = permute(itc, [3 1 2]);
                source_tf.msk(inside_idx, ...
                              1:length(source_tf.freq), :) = itc < maskthresh;
            case 'pow'
                source_tf.msk = ...
                    source_tf.avg.pow > maskthresh*max(source_tf.avg.pow); % threshold at some % of peak
            case 'momhigh'
                source_tf.msk(inside_idx,1,:) = mom' > maskthresh;
            case 'momlow'
                source_tf.msk(inside_idx,1,:) = mom' < maskthresh;
            case 'momabs'
                source_tf.msk(inside_idx,1,:) = abs(mom') > maskthresh;
            case 'momratio'
                for tt = 1:length(source_tf.time)
                    source_tf.msk(inside_idx,1,tt) = abs(mom(tt, :)) > ...
                        maskthresh*max(abs(mom(:)));
                end
            case 'pval'
                pval = reshape([source_tf.pval{:}], ...
                               length(source_tf.freq), ...
                               length(source_tf.time), ...
                               length(inside_idx));
                pval = permute(pval, [3 1 2]);
                source_tf.msk(inside_idx, ...
                              1:length(source_tf.freq),:) = pval < maskthresh;
            otherwise
                error('masktype unknown');
        end

        % plot it
        nmt_sourceplot(cfg,source_tf);

    else

        cfg = [];
        cfg.zlim = [-3 3];
        inside_idx = find(source_tf.inside);
        cfg.mripath = t1_template_fname;
        cfg.funparameter = 'stat';
        cfg.plottype = 'tf';  % ts if you want to view as time series instead
        cfg.atlas = ft_read_atlas(atlas_fname);
        nmt_sourceplot(cfg,ft_convert_units(source_tf,'mm'));

    end
end
