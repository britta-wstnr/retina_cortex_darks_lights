% Prepares data for Hilbert beamforming, will results in Fig. 4

plot_check = 1;

% Loop over subjects
for nn = 1:length(subjs)

    keep plot_check

    project_settings;
    display(sprintf('Processing subject %s', subjs{nn}));

    proc_dir      = fullfile(base_dir, subjs{nn});
    headmodel_dir = fullfile(base_dir, subjs{nn}, 'bem');
    mri_dir       = fullfile(base_dir, 'mri', ['s', subjs{nn}]);
    out_dir       = fullfile(base_dir, subjs{nn}, 'hilbert_out');

    if(exist(out_dir, 'dir') == 0)
        cd(fullfile(base_dir, subjs{nn}));
        mkdir hilbert_out
    end

    %% Volume

    load(fullfile(headmodel_dir, vol_fname));
    vol.basefile = [subjs{nn}, '_mags'];

    %% Load data

    load(fullfile(proc_dir, data_on_meg_fname));
    load(fullfile(proc_dir, data_off_meg_fname));
    load(fullfile(proc_dir, 'grad_coilacc.mat'));

    % select data
    cfg = [];
    cfg.channel = meg_channel_sel;
    cfg.latency = [-.5 .4];   % will be cut again - prevent filter artifacts
    dataon_clean  = ft_selectdata(cfg, dataon_clean);
    dataoff_clean = ft_selectdata(cfg, dataoff_clean);

    % downsample to 1000 Hz
    cfg = [];
    cfg.resamplefs = 1000;
    dataon_clean = ft_resampledata(cfg, dataon_clean);
    dataoff_clean = ft_resampledata(cfg, dataoff_clean);

    % fix coil accuracy
    dataon_clean.grad  = grad_coilacc;
    dataoff_clean.grad = grad_coilacc;

    %% DFT filter

    cfg = [];
    cfg.dftfilter = 'yes';
    cfg.dftfreq = [50 100 150];
    cfg.demean = 'yes';
    cfg.coilaccuracy = 2;
    dataon_clean = ft_preprocessing(cfg, dataon_clean);
    dataoff_clean = ft_preprocessing(cfg, dataoff_clean);

    %% Resliced MRI

    rmri_path = fullfile(mri_dir, ['s', subjs{nn}, '_resl.nii']);

    %% sourcemodel

    % template:
    load(source_template_fname);
    template_grid = ft_convert_units(sourcemodel, 'mm');

     % load individual source model
    load(fullfile(proc_dir, sourcem_fname));
    sourcemodel = ft_convert_units(sourcemodel, 'mm');

    %% load coregistration

    load(fullfile(proc_dir, coreg_fname));
    coreg = nuts.coreg;

    %% coregister data structure

    switch(dataon_clean.grad.type)
        case {'bti148','bti248'}
            origcoordsys = 'bti';
        case 'neuromag306'
            origcoordsys = 'neuromag';
        otherwise
            error('Fiducial convention unknown, add it to the cases above!')
    end

    % get coregistration tfm
    [coreg.mri2meg_tfm, coordsys] = ft_headcoordinates(...
        coreg.fiducials_mri_mm(3, :), ...
        coreg.fiducials_mri_mm(1, :), ...
        coreg.fiducials_mri_mm(2, :), origcoordsys);
    coreg.meg2mri_tfm = inv(coreg.mri2meg_tfm);

    % convert grad structure with nuts
    mri_grad = dataon_clean.grad;
    mri_grad = ft_convert_units(mri_grad,'mm');

    % transforms the grad_mm coordinates to mri mm coordinates
    mri_grad = ft_transform_sens(coreg.meg2mri_tfm, mri_grad);
    mri_grad.coordsys = 'spm';

    %% sanity check

    if(plot_check)
        figure; hold on     % plot all objects in one figure
        ft_plot_vol(vol, 'edgecolor', 'none', 'facecolor', 'blue')
        alpha 0.1          % make the surface transparent
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
        a=ft_plot_sens(mri_grad, 'style', 'r*');
        set(a, 'linewidth', 2)
        view([90 0 0])
    end

    %% make copy of data to prevent coregistering already coregistered data

    dataon_coreg = dataon_clean; dataoff_coreg = dataoff_clean;
    dataon_coreg.grad = mri_grad;
    dataoff_coreg.grad = mri_grad;

    %% LOAD LEADFIELD

    load(fullfile(headmodel_dir, lf_name));

    %% Hilberting ON and OFF

    display('preprocessing done, starting Hilberting')

    stim_def = 'on';
    display('Starting ON Analysis')
    hilbert_meg

    display('Starting OFF Analysis')
    stim_def = 'off';
    hilbert_meg
end
