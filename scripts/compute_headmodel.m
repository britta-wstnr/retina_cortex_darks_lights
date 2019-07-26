% Compute BEM headmodels from MRI

project_settings;

% seg method that works for all is 'fieldtrip'
seg_method = 'fieldtrip' % 'load_spm', 'fieldtrip' or 'spm_newseg';
check_mesh = true;

nn = 1;  % doing this manually, as looping is not feasible with checking every
         % step
subj_id = subjs{nn};

subj_dir      = fullfile(base_dir, subjs{nn});
headmodel_dir = fullfile(base_dir, subjs{nn}, 'bem');
mri_dir       = fullfile(base_dir, 'mri', ['s', subjs{nn}]);
mri_file_resl = fullfile(mri_dir, ['s', subjs{nn}, '_resl.nii']);

openmeeg_filename = [subj_id];

[mri_file, update_resl, resl_new] = get_mri_file(subj_id, mri_dir)
if update_resl
    mri_file_resl = resl_new;
end

%% Create directories if necessary

if(exist(headmodel_dir, 'dir')==0)
    cd(fullfile(base_dir, subj_id));
    mkdir bem
end

%% Reslicing

cd(mri_dir)
if(exist(mri_file_resl, 'file')==0)
    reslice_nii(mri_file, mri_file_resl, 1);
end

mri = ft_read_mri(mri_file_resl);
mri.coordsys = 'spm';

ft_sourceplot([], mri);

%% spm8 newseg or fieldtrip segmentation

switch seg_method

    case 'spm_newseg'
        seg = create_bem_segmentation('t1filename', mri_file_resl);
        seg.transform = mri.transform;
        seg.dim = mri.dim;
        seg.unit = mri.unit;

        filename = fullfile(headmodel_dir, 'spm_seg.mat')
        save(filename, 'seg');

    case 'fieldtrip'
        cfg = [];
        cfg.output = {'gray','white','csf','skull','scalp'};
        seg = ft_volumesegment(cfg,mri);

        filename = fullfile(headmodel_dir, 'ft_seg_m.mat')
        save(filename, 'seg');

    case 'load_spm'

        load(fullfile(headmodel_dir, 'spm_seg.mat'));
end

%% Check segmentation

% add anatomical information to the segmentation
seg.transform = mri.transform;
seg.anatomy   = mri.anatomy;

% plot it:
cfg = [];
switch seg_method
    case {'spm_newseg' 'load_spm'}
        cfg.funparameter = 'brain'
    case 'fieldtrip'
        cfg.funparameter = 'gray';
end
ft_sourceplot(cfg,seg);

cfg.funparameter = 'skull';
ft_sourceplot(cfg,seg);

cfg.funparameter = 'scalp';
ft_sourceplot(cfg,seg);

%% Prepare mesh

cfg = [];
cfg.tissue = {'scalp', 'skull', 'brain'};
cfg.method = 'iso2mesh';    % 'projectmesh';
cfg.numvertices = 10000;    % We'll decimate later
bnd = ft_prepare_mesh(cfg, seg);

% Decimate
[bnd(1).pos, bnd(1).tri] = meshresample(bnd(1).pos, ...
                                        bnd(1).tri, ...
                                        1000/size(bnd(1).pos, 1));   % scalp
[bnd(2).pos, bnd(2).tri] = meshresample(bnd(2).pos, ...
                                        bnd(2).tri, ...
                                        2000/size(bnd(2).pos, 1));   % skull
[bnd(3).pos, bnd(3).tri] = meshresample(bnd(3).pos, ...
                                        bnd(3).tri, ...
                                        3000/size(bnd(3).pos, 1));   % brain

% Check and repair individual meshes using iso2mesh
for ii = 1:length(bnd)
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, ...
                                                 bnd(ii).tri, 'dup');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, ...
                                                 bnd(ii).tri, 'isolated');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, ...
                                                 bnd(ii).tri, 'deep');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, ...
                                                 bnd(ii).tri, 'meshfix');
end

% Ensure no overlaps
bnd = decouplesurf(bnd);

%% plot the output as sanity check

switch seg_method

    case 'spm_newseg'
        brain = bnd(3);

        % smooth brain surface
        [conn, connnum, count] = meshconn(brain(1).tri, size(brain(1).pos, 1));
        brain(1).pos = smoothsurf(brain(1).pos, [], conn, 10, 0.9, 'lowpass');

        display(subj_id)   % to keep track while looking at the plots

        if check_mesh == true
            figure;
            ft_plot_mesh(brain(1), 'facecolor', 'red'); hold on; pause
            ft_plot_mesh(bnd(2), 'facealpha', 0.25, 'facecolor', 'blue');
            hold on; pause             % skull
            ft_plot_mesh(bnd(1), 'facealpha', 0.25);% pause
            bnd(3) = brain(1);
        end

    case 'fieldtrip'
        if check_mesh == true
            figure;
            ft_plot_mesh(bnd(3), 'facecolor', 'red'); hold on; pause
            ft_plot_mesh(bnd(2), 'facealpha', 0.25, 'facecolor', 'blue');
            hold on; pause             % skull
            ft_plot_mesh(bnd(1), 'facealpha', 0.25);
        end
end

if(0)
    % smooth scalp surface if needed
    scalp = bnd(1);
    [conn,connnum,count] = meshconn(scalp(1).tri, size(scalp(1).pos,1));
    scalp(1).pos=smoothsurf(scalp(1).pos, [], conn, 10, 0.9, 'lowpass');

    figure;
    ft_plot_mesh(bnd(3), 'facecolor', 'red'); hold on; pause
    ft_plot_mesh(bnd(2), 'facealpha', 0.25, 'facecolor', 'blue');
    hold on; pause             % skull
    ft_plot_mesh(scalp, 'facealpha', 0.25);

    figure;
    ft_plot_mesh(scalp, 'facecolor', 'blue');
    figure;
    ft_plot_mesh(bnd(1), 'facecolor', 'green');

    bnd(1) = scalp(1);
end


%% Create the BEM volume and save it for later

vol.bnd = bnd;
vol.cond = [0.33 0.0041 0.33]; % SI units, ignore CSF
vol.type = 'openmeeg';
vol.basefile = [openmeeg_filename, '_3layer'];
% following files in here can be reused: hm.bin, hm_inv.bin, dsm.bin
vol.path = [headmodel_dir, '/openmeeg_out'];
vol = ft_convert_units(vol,'mm');    % Convert bnd to SI units

vol.segmentationmethod = seg_method;

filename = fullfile(headmodel_dir, vol_fname);
save(filename, 'vol');
