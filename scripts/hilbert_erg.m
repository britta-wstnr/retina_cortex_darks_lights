% Compute Hilbert transform of ERG data

%%  script-speficic specifications
stim_def = 'on';  % this is gonna be used by nemo_hilbert_sensorlevel
baseline = 'normal'

for nn = 1:length(subjs)

    keep fieldtrip_path fieldtrip_path_nmt subjs stim_def nn baseline
    project_settings;

    proc_dir      = fullfile(base_dir, subjs{nn});
    headmodel_dir = fullfile(base_dir, subjs{nn}, 'bem');
    mri_dir       = fullfile(base_dir, 'mri', ['s', subjs{nn}]);
    out_dir       = fullfile(base_dir, subjs{nn}, 'hilbert_out');

    if(exist(out_dir, 'dir')==0)
        cd(fullfile(base_dir, subjs{nn}));
        mkdir hilbert_out
    end

    %% Load data
    load(fullfile(proc_dir, data_on_erg_fname));
    load(fullfile(proc_dir, data_off_erg_fname));

    cfg = [];
    cfg.channel = erg_channel_sel;
    cfg.latency = [-.5 .4];   % will be cut again - prevent filter artifacts
    dataon_clean  = ft_selectdata(cfg, dataon_clean2);
    dataoff_clean = ft_selectdata(cfg, dataoff_clean2);

    % downsample to 1000 Hz
    cfg = [];
    cfg.resamplefs = 1000;
    dataon_clean = ft_resampledata(cfg, dataon_clean);
    dataoff_clean = ft_resampledata(cfg, dataoff_clean);


    %% DFT filter

    cfg = [];
    cfg.dftfilter = 'yes';
    cfg.dftfreq = [50 100 150];
    cfg.demean = 'yes';
    dataon_clean = ft_preprocessing(cfg, dataon_clean);
    dataoff_clean = ft_preprocessing(cfg, dataoff_clean);

    %% Run for ON and OFF

    stim_def = 'on';   % this is gonna be used by nemo_hilbert_sensorlevel
    nemo_hilbert_sensorlevel

    stim_def = 'off';   % this is gonna be used by nemo_hilbert_sensorlevel
    nemo_hilbert_sensorlevel

end