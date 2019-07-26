% Compute the leadfields for source reconstruction

project_settings;

plot_check = 1;

for nn = 1:length(subjs)

    proc_dir      = fullfile(base_dir, subjs{nn});
    headmodel_dir = fullfile(base_dir, subjs{nn}, 'bem');
    mri_dir       = fullfile(base_dir, 'mri', ['s', subjs{nn}]);

    %% Load volume

    load(fullfile(headmodel_dir, vol_fname);
    vol.basefile = [subjs{nn}, '_mags'];

    %% Data

    load(fullfile(proc_dir, data_on_meg_fname);

    %% MRI s0001

    rmripath = fullfile(mri_dir, ['s', subjs{nn}, '_resl.nii']);
    reslice_nii(mripath, rmripath)

    clear mripath   % to be sure not to work with the wrong MRI

    %% sourcemodel

    % template:
    load(source_template_fname);
    template_grid = ft_convert_units(sourcemodel, 'mm');

    % warp if not done yet:
    if(0)
        mri = ft_read_mri(rmripath);
        mri.coordsys = 'spm';

        cfg                = [];
        cfg.grid.warpmni   = 'yes';
        cfg.grid.template  = template_grid;
        cfg.grid.nonlinear = 'yes'; % use non-linear normalization
        cfg.mri            = mri;
        sourcemodel = ft_prepare_sourcemodel(cfg);
        sourcemodel = ft_convert_units(sourcemodel, 'mm');

        save(fullfile(proc_dir, sourcem_fname), 'sourcemodel');
    else
        load(fullfile(proc_dir, sourcem_fname));
        sourcemodel = ft_convert_units(sourcemodel, 'mm');
    end

    %% do manual coregistration of MRI in nutmeg

     if(exist(fullfile(proc_dir, 'nuts_nii.mat'), 'file')==0)

        cd (fullfile(mri_dir));
        nm
        pause
        cd(proc_dir)
        global nuts
        nuts.fig = 1  % erase the figure so MATLAB doesn't freak when reloading
        save nuts_nii nuts

    else
        load(fullfile(proc_dir, coreg_fname));
    end

    coreg = nuts.coreg;

    %% get coregistration tfm
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

    dataon_coreg = dataon_clean;
    dataon_coreg.grad = mri_grad;

    %% COMPUTE LEADFIELD

    cfg             = [];
    cfg.grid        = sourcemodel;
    cfg.reducerank  = 'no';
    cfg.normalize   = 'no';  % no LF norm but normalizing BF weights instead
    cfg.headmodel   = vol;
    cfg.channel     = meg_channel_sel;
    cfg.grad        = dataon_coreg.grad;
    leadfield       = ft_prepare_leadfield(cfg);

    save(fullfile(headmodel_dir, lf_fname), 'leadfield');

end