% Compute Hilbert amplitude and ITC on sensor-level data, wrapper script

% this is a version of
% https://github.com/meeg-cfin/nemolab/blob/master/hilbertsourceloc/nemo_hilbert_sensorlevel.m

% variable space needs to have loaded:
% data (renamed below!)

% output of this script can either be plotted with nemo_plot_hilbertchannel
% (make sure to specify the frequency bands as lower and upper limits per
% band!) or with fieldtrip's own plotting functions, ft_singleplotTFR and
% ft_multiplotTFR (if topo-data) - then make sure to specify data.freq
% containing the mean of each frequency band.

% specs for firws filter:
tswidth = 6;

baselinewindow = [-.15 -.05];
toilim  = [-.15 .15]

stim = stim_def ; %'on';
filtertype = 'firws'; % 'but' or 'firws'

if(strcmp(stim, 'on'))
    data = dataon_clean;  % data already loaded - adjust for your needs.
                          % This is preprocessed + downsampled sensorlevel data.
    display('Analyzing ON data');
elseif(strcmp(stim, 'off'))
    data = dataoff_clean;
    display('Analyzing OFF data');
else
    error('unknown stimulus condition!')
end

saving = 1;     % should data be saved?
plotting = 0;   % should results be plotted?

%% run through filter bank, and obtain timelock and hilbert transform

% preallocation
timelockbp = cell(1,size(freq_bands, 1));
data_hilb_bp  = cell(1,size(freq_bands, 1));
filterdefinitions = cell(1,size(freq_bands, 1));

% using parfor here has some weird effects on filtering feedback etc
for ii = 1:size(freq_bands,1)
    cfg=[];
    cfg.demean = 'yes';
    cfg.baselinewindow = baselinewindow;
    if(strcmp(filtertype, 'but'))
        cfg.bpfilter = 'yes'; % NB: default butterworth for quick testing;
                              % specify more advanced filter for real analysis!
        cfg.bpfreq = freq_bands(ii,:);
    elseif(strcmp(filtertype, 'firws'))
        cfg.hpfilter = 'yes';
        cfg.hpfreq = freq_bands(ii,1);
        cfg.hpfilttype = filtertype;
        cfg.hpfiltdf = tswidth;
        cfg.hpfiltdir =  'onepass-zerophase';
        cfg.hpfiltwintype = 'hann';
        cfg.plotfiltresp = 'no';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = freq_bands(ii,2);
        cfg.lpfilttype = filtertype;
        cfg.lpfiltdf = tswidth;
        cfg.lpfiltdir = 'onepass-zerophase';
        cfg.lpfiltwintype = 'hann';
    else
        error('I do not know the specified filter type, please add!')
    end
    databp = ft_preprocessing(cfg,data);

    % filtering first then snipping to shorter time interval avoids edge
    % artifacts
    cfg = [];
    cfg.toilim = toilim;
    databp = ft_redefinetrial(cfg,databp);

    cfghilb=[];
    cfghilb.hilbert = 'complex';
    data_hilb_bp{ii} = ft_preprocessing(cfghilb,databp);

end
clear databp


%% apply this now to the complex-valued hilbert transform of the data

saveRAM = 0;
tfstats = 1;

disp(['Hilbert-analysis started at ' datestr(now)])
tstart = tic;
for ii = 1:size(freq_bands, 1)
    % perform stats on sensorlevel
    nemo_hilbertstats_sensorlevel
end
disp(['Hilbert-analysis finished at ' datestr(now) ' and took ' ...
      num2str(toc(tstart)/60) ' minutes to run']);

if(saveRAM)
    clear data_hilb_bp
end

%% assembles source_hilb{:} into a composite source_tf structure

data_tf = data_hilbert{1};
data_tf.avg = zeros(numel(data_hilbert{1}.label), size(freq_bands,1), ...
                    size(data_hilbert{1}.time, 2));
data_tf.itc = data_tf.avg;
data_tf.stat = data_tf.avg;
data_tf.pval = data_tf.avg;

data_tf.freq_bands = freq_bands;

% This is needed if output should be plotted with fieldtrip's plotting
% functions - ft_singleplotTFR or ft_multiplotTFR
data_tf.freq = mean(data_tf.freq_bands');

for jj = 1:size(freq_bands,1)

    data_tf.avg(:, jj, :) = reshape(data_hilbert{jj}.avg, ...
                                    numel(data_hilbert{1}.label), 1, ...
                                    size(data_hilbert{1}.time, 2));
    data_tf.itc(:, jj, :) = reshape(data_hilbert{jj}.itc, ...
                                    numel(data_hilbert{1}.label), 1, ...
                                    size(data_hilbert{1}.time, 2));

    if(tfstats)
        data_tf.stat(:, jj, :) = reshape(data_hilbert{jj}.stat, ...
                                         numel(data_hilbert{1}.label), 1, ...
                                         size(data_hilbert{1}.time, 2));
        data_tf.pval(:, jj, :) = reshape(data_hilbert{jj}.pval, ...
                                         numel(data_hilbert{1}.label), 1, ...
                                         size(data_hilbert{1}.time, 2));
    end
end

data_tf.dimord = 'chan_freq_time';
if(saving)
    if(strcmp(stim, 'on'))
        save(fullfile(out_dir, 'erg_ON_rightstats.mat'), 'data_tf');
        display('Saving ON data');
    elseif(strcmp(stim, 'off'))
        save(fullfile(out_dir, 'erg_OFF_rightstats.mat'), 'data_tf');
        display('Saving OFF data');
    end
end
